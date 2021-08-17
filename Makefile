include config.mk

OBJFILES	=	obj/IHSOptimizer.o obj/Main.o obj/CPLEXInterface.o obj/Timer.o obj/ArgsParser.o obj/GlobalConfig.o obj/Util.o

CPPFLAGS	+=	-O3 -std=c++14 -pthread -fPIC -m64 \
-D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS -MMD

WARNS	=	-pedantic -Wall -Wextra -Wno-ignored-attributes -Wno-reorder -Wno-sign-compare

CPPFLAGS	+= -fno-strict-aliasing -fexceptions -DIL_STD -Ioptions
IP_LNDIRS	=	-L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
IP_LNFLAGS	=	-lconcert -lilocplex -lcplex -lm -ldl -lpthread
IP_INCLUDES	=	-I$(CPLEXDIR)/include -I$(CONCERTDIR)/include

WASP_INCLUDES	=	-I./wasp/src


EXE	=	bin/asp-hs


.PHONY: all clean wasp dirs
all: options/regenerate $(EXE)

options/%.cpp:
	@rm -f options/regenerate

options/%.h:
	@rm -f options/regenerate

options/regenerate: options/params.csv options/generate.py
	@echo "Generating options"
	@cd options && ./generate.py
	@touch options/regenerate

clean:
	rm -rf obj bin
	make -s -C ./wasp clean

dirs:
	@-mkdir -p obj
	@-mkdir -p bin

wasp:
	make -s -C ./wasp BUILD=release

$(EXE):	dirs $(OBJFILES) wasp
	@echo "-> linking $@"
	@g++ $(CPPFLAGS) $(OBJFILES) $(shell find ./wasp/build/release/ -name *.o | grep -v main) $(IP_LNDIRS) $(IP_LNFLAGS) -o $@

obj/%.o: src/%.cpp
	@echo "-> compiling $@"
	@g++ $(CPPFLAGS) $(WASP_INCLUDES) $(IP_INCLUDES) $(WARNS) -c $< -o $@

-include $(shell [ -e bin ] && find obj -type f -name "*.d")

########## Tests

TESTS_DIR = wasp/tests

TESTS_TESTER = $(TESTS_DIR)/pyregtest.py

BINARY_ASPHS = $(EXE)
TESTS_COMMAND_waspweak = $(BINARY_ASPHS) - --silent --disjoint-hs --reducedcosts --clarke-lp-lb --ip-seeding

TESTS_COMMAND_WeakConstraints = $(TESTS_COMMAND_waspweak)

TESTS_CHECKER_WeakConstraints = $(TESTS_DIR)/weakConstraints.checker.py

TESTS_REPORT_text = $(TESTS_DIR)/text.report.py

TESTS_DIR_asp_WeakConstraints = $(TESTS_DIR)/asp/weakConstraints
TESTS_SRC_asp_WeakConstraints = $(sort $(shell find $(TESTS_DIR_asp_WeakConstraints) -name '*.test.py'))
TESTS_OUT_asp_WeakConstraints = $(patsubst %.test.py,%.test.py.text, $(TESTS_SRC_asp_WeakConstraints))

test: wasp/tests/asp/weakConstraints

wasp/tests/asp/weakConstraints: $(TESTS_OUT_asp_WeakConstraints)

$(TESTS_OUT_asp_WeakConstraints):
	@$(TESTS_TESTER) "$(TESTS_COMMAND_WeakConstraints)" $(patsubst %.test.py.text,%.test.py , $@) $(TESTS_CHECKER_WeakConstraints) $(TESTS_REPORT_text)
