#! /usr/bin/env python
"""
multi_dict.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from utils.error import MyError
from utils.messages import TestMsg, ErrMsg

def SetupMultiDict(top_dict, keys, item, list_end=True): #{
  AddToMultiDict(top_dict, keys, item, list_end)
#} end def

def AddToMultiDict(top_dict, keys, item, list_end=True, doextend=False): #{
  curr_dict = top_dict
  for key in keys[:-1]: #{
    curr_dict = AddToMultiDictLevel(curr_dict, key)
  #} end for
  if (list_end): #{
    bottom_list = AddToMultiDictLevel(curr_dict, keys[-1], bottom=True)
    if (doextend and isinstance(item, list)): #{
      bottom_list.extend(item)
    else:
      bottom_list.append(item)
    #} end if
    #print ("Appending %s\n  Keys: %s\n  Now: %s" %
    #  (item, ",".join(keys), ",".join(bottom_list)))
  elif (keys[-1] in curr_dict): #{
    if (curr_dict[keys[-1]] != item): #{
      raise MultiDictError("value for keys %s already exists: %s != %s" %
        (",".join(map(str, keys)), curr_dict[keys[-1]], item))
    #} end if
  else:
    curr_dict[keys[-1]] = item
  #} end if
  #chromA_dict = SetupMultiDict(self.barnacle_results, chromA)
  #chromB_dict = SetupMultiDict(chromA_dict, chromB)
  #posA_dict = SetupMultiDict(chromB_dict, posA)
  #posB_list = SetupMultiDict(posA_dict, posB, bottom=True)
  #posB_list.append(event)
#} end def

def AddToMultiDictLevel(curr_dict, new_key, bottom=False): #{
  if (new_key not in curr_dict): #{
    if (bottom): #{
      curr_dict[new_key] = list()
    else:
      curr_dict[new_key] = dict()
    #} end if
  #} end if
  return curr_dict[new_key]
#} end def

def CheckMultiDict(top_dict, keys): #{
  curr_dict = top_dict
  for key in keys:
    if (key not in curr_dict): #{
      return False
    #} end if
    curr_dict = curr_dict[key]
  #} end for
  return True
#} end def

def AccessMultiDict(top_dict, keys): #{
  curr_dict = top_dict
  for key in keys[:-1]: #{
    curr_dict = curr_dict[key]
  #} end for
  return curr_dict[keys[-1]]
#} end def

def TraverseMultiDict(curr_dict, sort=False, sort_key=None): #{
  if (isinstance(curr_dict, dict)): #{
    if (sort): #{
      key_list = sorted(curr_dict.keys(), key = sort_key)
      if (None != sort_key and hasattr(sort_key, "NextLevel")): #{
        sort_key.NextLevel()
      #} end if
    else:
      key_list = curr_dict.keys()
    #} end if
    for key in key_list: #{
      for elem in TraverseMultiDict(curr_dict[key], sort,
          sort_key=sort_key): #{
        yield elem
      #} end for
    #} end for
    if (sort and None != sort_key and hasattr(sort_key, "PrevLevel")): #{
      sort_key.PrevLevel()
    #} end if
  elif (isinstance(curr_dict, list)):
    if (sort): #{
      curr_list = sorted(curr_dict, key=sort_key)
    else:
      curr_list = curr_dict
    #} end if
    for elem in curr_list: #{
      yield elem
    #} end for
  else:
    yield curr_dict
  #} end if
#} end def

def StringifyMultiDict(curr_dict, sort=False, sort_key=None): #{
  if (isinstance(curr_dict, str)): #{
    return curr_dict
  #} end if
  if (isinstance(curr_dict, list)): #{
    if (sort): #{
      return ",".join(sorted(curr_dict, key=sort_key))
    else:
      return ",".join(curr_dict)
    #} end if
  #} end if
  if (sort): #{
    key_list = sorted(curr_dict.keys(), key=sort_key)
    if (None != sort_key and hasattr(sort_key, "NextLevel")): #{
      sort_key.NextLevel()
    #} end if
  else:
    key_list = curr_dict.keys()
  #} end if
  output_list = list()
  for key in key_list: #{
    output_list.append("%s(%s)" % (key,
      StringifyMultiDict(curr_dict[key], sort=sort, sort_key=sort_key)))
  #} end for
  if (sort and None != sort_key and hasattr(sort_key, "PrevLevel")): #{
    sort_key.PrevLevel()
  #} end if
  return ",".join(output_list)
#} end def

def RemoveFromMultiDict(top_dict, keys, clear_empty=True): #{
  try:
    RemoveFromMultiDictLevel(top_dict, keys, clear_empty=clear_empty)
  except MultiDictError, e:
    raise MultiDictError("%s. Keys: %s" % (e.msg, ",".join(keys)))
  #} end try
#} end def

def RemoveFromMultiDictLevel(curr_dict, keys, clear_empty=True): #{
  if (1 == len(keys)): #{
    if (isinstance(curr_dict, dict) or isinstance(curr_dict, list)): #{
      del curr_dict[keys[0]]
    else:
      raise MultiDictError("Error removing item: too many keys")
    #} end if
  else:
    RemoveFromMultiDict(curr_dict[keys[0]], keys[1:], clear_empty=clear_empty)
    if ((isinstance(curr_dict, dict) or isinstance(curr_dict, list)) and
        clear_empty and 0 == len(curr_dict)): #{
      del curr_dict[keys[0]]
    #} end if
  #} end if
#} end def

class MultiDictCls: #{
  def __init__(self, input_str, depth, list_end=True): #{
    #ErrMsg("  Initializing multi-dict. Input: %s" % input_str)
    self.top_dict = dict()
    self.input_str = input_str
    self.depth = depth
    self.list_end = list_end
    if ("none" != input_str.lower()): #{
      if (not input_str.endswith(")")): #{
        raise MultiDictError("Misformed multi-dict input string, must "
          "end with a closing parenthesis: %s" % input_str)
      #} end if
      self.ParseInputString()
    #} end if
  #} end def

  def __getitem__(self, key): #{
    return self.top_dict[key]
  #} end def

  def __iter__(self): #{
    return self.top_dict.iterkeys()
  #} end def

  def __len__(self): #{
    return len(self.top_dict)
  #} end def

  def __str__(self): #{
    return self.OutputString()
  #} end def

  def ParseInputString(self): #{
    self.Initialize()
    prev_pos = -1
    for delim_pos in self.delim_positions: #{
      curr_value = self.input_str[prev_pos+1:delim_pos]
      curr_delim = self.input_str[delim_pos]
      TestMsg("Curr: \"%s\" \"%s\" (%i)" % (curr_value, curr_delim,
        self.curr_depth))
      # an opening parenthesis means going deeper
      if ("(" == curr_delim): #[
        self.IncreaseDepth(curr_value, delim_pos)
      # a comma means going on to the next item
      elif ("," == curr_delim):
        self.MaintainDepth(curr_value, delim_pos)
      # a closing parenthesis means getting less deep
      elif (")" == curr_delim): #{
        self.DecreaseDepth(curr_value, delim_pos)
      else:
        raise MultiDictError("Unrecognized delimiter \"%s\" "
          "encountered while parsing multi-dict string: %s" %
          (curr_delim, self.input_str))
      #} end if
      prev_pos = delim_pos
    #} end for
  #} end def

  def Initialize(self): #{
    # clear the dictionary
    self.top_dict = dict()
    # initialize the keys
    self.keys = [None for i in range(self.depth-1)]
    # get the positions of all parentheses and commas
    self.delim_positions = [i for i,c in
      enumerate(self.input_str) if c in "(,)"]
    TestMsg("DELIMS: %s" % ",".join(map(str, self.delim_positions)))
    # start at depth zero
    self.curr_depth = 0
  #} end def

  def IncreaseDepth(self, curr_value, delim_pos): #{
    TestMsg("INCREASING DEPTH!")
    self.keys[self.curr_depth] = curr_value
    self.curr_depth += 1
    # make sure that we have not gone too deep
    if ((self.depth-1) < self.curr_depth): #{
      raise MultiDictError("Misplaced opening parenthesis "
        "encountered while parsing multi-dict string: %s (%i)" %
        (self.input_str, delim_pos))
    #} end if
  #} end def

  def MaintainDepth(self, curr_value, delim_pos): #{
    TestMsg("MAINTAINING DEPTH!")
    if ((self.depth-1) == self.curr_depth): #{
      SetupMultiDict(self.top_dict, self.keys, curr_value, self.list_end)
    # if we are not at the max depth, this should only occur after
    # a closing parentheis
    elif (")" != self.input_str[delim_pos-1]):
      raise MultiDictError("Misplaced comma encountered while "
        "parsing multi-dict string: %s (%i)" % (self.input_str, delim_pos))
    #} end if
  #} end def

  def DecreaseDepth(self, curr_value, delim_pos): #{
    TestMsg("DECREASING DEPTH!")
    # if we are at the max depth, add the current value to the dictionary
    if ((self.depth-1) == self.curr_depth): #{
      SetupMultiDict(self.top_dict, self.keys, curr_value, self.list_end)
    # if we are not at the max depth, this should only occur after
    # another closing parentheis
    elif (")" != self.input_str[delim_pos-1]):
      raise MultiDictError("Misplaced closing parenthesis "
        "encountered while parsing multi-dict string: %s (%i)" %
        (self.input_str, delim_pos))
    #} end if
    self.curr_depth -= 1
    # make sure that we have not gotten too shallow
    if (0 > self.curr_depth): #{
      raise MultiDictError("Misplaced closing parenthesis "
        "encountered while parsing multi-dict string: %s (%i)" %
        (self.input_str, delim_pos))
    #} end def
  #} end def

  def Add(self, keys, item): #{
    AddToMultiDict(self.top_dict, keys, item, self.list_end)
  #} end def

  def OutputString(self, sort=False, sort_key=None): #{
    return StringifyMultiDict(self.top_dict, sort=sort, sort_key=sort_key)
  #} end def
#} end class

class MultiDictSortKeyCls: #{
  def __init__(self): #{
    self.level = 0
  #} end def

  def NextLevel(self): #{
    self.level += 1
  #} end def

  def PrevLevel(self): #{
    self.level -= 1
  #} end def

  def __call__(self, arg1): #{
    return arg1
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class MultiDictError(MyError): #{
  pass
#} end class
