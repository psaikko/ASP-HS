
# parsetab.py
# This file is automatically generated. Do not edit.
_tabversion = '3.2'

_lr_method = 'LALR'

_lr_signature = '7^\x90}@\xb1\xb2\xbeD\x8f\xbe3i\x9f\xcbS'
    
_lr_action_items = {'COUNT':([27,],[34,]),'PAR_CLOSE':([21,22,23,24,25,26,39,],[-24,32,-26,-27,-28,-29,-25,]),'GREAT':([4,5,7,9,32,63,],[-19,-23,-21,-22,-20,67,]),'CURLY_OPEN':([34,35,36,37,],[40,41,42,43,]),'MIN':([27,],[37,]),'MAX':([27,],[35,]),'SUM':([27,],[36,]),'LESS':([40,41,42,43,51,],[45,45,45,45,45,]),'CURLY_CLOSE':([44,46,47,48,49,57,67,],[50,-16,53,54,55,-17,-18,]),'COLON':([52,],[58,]),'ID':([0,1,3,8,10,11,12,13,18,28,29,31,33,58,],[-1,7,-2,7,7,7,-5,23,7,7,-4,-3,23,7,]),'LEQ':([15,50,53,54,55,],[27,56,59,60,61,]),'NUM':([8,10,13,28,33,45,56,59,60,61,],[15,15,24,15,24,52,62,64,65,66,]),'PAR_OPEN':([4,5,7,9,],[13,-23,-21,-22,]),'STR':([0,1,3,8,10,11,12,13,18,28,29,31,33,58,],[-1,5,-2,5,5,5,-5,25,5,5,-4,-3,25,5,]),'NOT':([8,10,28,],[18,18,18,]),'COMMA':([4,5,7,9,14,16,17,19,21,22,23,24,25,26,30,32,38,39,44,46,47,48,49,57,62,64,65,66,67,],[-19,-23,-21,-22,-6,28,-10,28,-24,33,-26,-27,-28,-29,-11,-20,-7,-25,51,-16,51,51,51,-17,-12,-14,-13,-15,-18,]),'DERIVES':([0,1,2,3,4,5,6,7,9,12,20,29,31,32,],[-1,8,10,-2,-19,-23,-8,-21,-22,-5,-9,-4,-3,-20,]),'OR':([0,1,2,3,4,5,6,7,8,9,10,11,12,13,18,20,28,29,31,32,33,58,],[-1,9,11,-2,-19,-23,-8,-21,9,-22,9,9,-5,26,9,-9,9,-4,-3,-20,26,9,]),'DOT':([2,4,5,6,7,9,14,16,17,19,20,30,32,38,62,64,65,66,],[12,-19,-23,-8,-21,-22,-6,29,-10,31,-9,-11,-20,-7,-12,-14,-13,-15,]),'$end':([0,1,3,12,29,31,],[-1,0,-2,-5,-4,-3,]),}

_lr_action = { }
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = { }
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'body':([8,10,],[16,19,]),'head':([1,],[2,]),'terms':([13,],[22,]),'set_terms':([40,41,42,43,],[44,47,48,49,]),'term':([13,33,],[21,39,]),'rule':([1,],[3,]),'set_term':([40,41,42,43,51,],[46,46,46,46,57,]),'literal':([8,10,28,],[14,14,38,]),'program':([0,],[1,]),'atom':([1,8,10,11,18,28,58,],[6,17,17,20,30,17,63,]),'id':([1,8,10,11,18,28,58,],[4,4,4,4,4,4,4,]),}

_lr_goto = { }
for _k, _v in _lr_goto_items.items():
   for _x,_y in zip(_v[0],_v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = { }
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> program","S'",1,None,None,None),
  ('program -> <empty>','program',0,'p_program_emptyline','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',11),
  ('program -> program rule','program',2,'p_program_rec','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',15),
  ('rule -> head DERIVES body DOT','rule',4,'p_rule','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',19),
  ('rule -> DERIVES body DOT','rule',3,'p_rule_constraint','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',23),
  ('rule -> head DOT','rule',2,'p_rule_facts','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',27),
  ('body -> literal','body',1,'p_body','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',31),
  ('body -> body COMMA literal','body',3,'p_body_rec','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',35),
  ('head -> atom','head',1,'p_head','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',39),
  ('head -> head OR atom','head',3,'p_head_rec','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',43),
  ('literal -> atom','literal',1,'p_literal','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',47),
  ('literal -> NOT atom','literal',2,'p_literal_neg','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',51),
  ('literal -> NUM LEQ COUNT CURLY_OPEN set_terms CURLY_CLOSE LEQ NUM','literal',8,'p_literal_count','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',55),
  ('literal -> NUM LEQ SUM CURLY_OPEN set_terms CURLY_CLOSE LEQ NUM','literal',8,'p_literal_sum','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',66),
  ('literal -> NUM LEQ MAX CURLY_OPEN set_terms CURLY_CLOSE LEQ NUM','literal',8,'p_literal_max','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',77),
  ('literal -> NUM LEQ MIN CURLY_OPEN set_terms CURLY_CLOSE LEQ NUM','literal',8,'p_literal_min','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',88),
  ('set_terms -> set_term','set_terms',1,'p_set_terms','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',99),
  ('set_terms -> set_terms COMMA set_term','set_terms',3,'p_set_terms_rec','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',103),
  ('set_term -> LESS NUM COLON atom GREAT','set_term',5,'p_set_term','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',108),
  ('atom -> id','atom',1,'p_atom_prop','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',112),
  ('atom -> id PAR_OPEN terms PAR_CLOSE','atom',4,'p_atom','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',116),
  ('id -> ID','id',1,'p_id_ID','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',120),
  ('id -> OR','id',1,'p_id_OR','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',124),
  ('id -> STR','id',1,'p_id_STR','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',128),
  ('terms -> term','terms',1,'p_terms','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',132),
  ('terms -> terms COMMA term','terms',3,'p_terms_rec','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',136),
  ('term -> ID','term',1,'p_term','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',140),
  ('term -> NUM','term',1,'p_term_num','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',144),
  ('term -> STR','term',1,'p_term_string','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',148),
  ('term -> OR','term',1,'p_term_OR','/home/marco/workspace/wasp/./black_box_test/../RewriteAggregates/rewrite_aggregates.py',152),
]
