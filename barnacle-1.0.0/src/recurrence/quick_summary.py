#! /usr/bin/env python
"""
quick_summary.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, re, sys, time, traceback

# import custom modules
from version import VERSION
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand, PowerSet)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  FileBoxCls)
from utils.library_iterator import LibIteratorCls
from parsers.candidate_group_parser import CandidateGroupParserCls

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"

EVENT_TYPES = [
  ("standard fusion", "fus"),
  ("mixed-sense fusion", "fus.mix_dirs"),
  ("total fusion", "fus.all"),
  ("ptd", "ptd"),
  ("itd-2f", "itd"),
  ("itd-2f multi", "itd.multi_exon"),
  ("itd-1f", "itd.edge_gap"),
  ("itd-1f multi", "itd.edge_gap.multi_exon"),
  ("itd-0f", "itd.full_ctg_dup"),
  ("itd-0f multi", "itd.full_ctg_dup.multi_exon"),
  ("itd-flanked", "itd.flanked"),
  ("itd-all", "itd.all"),
]

class PredictionSummarizerCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    #self.event_counters = dict(
    #  [(event_type[0], EventCounterCls(event_type[0], event_type[1])) for
    #    event_type in EVENT_TYPES])
    self.event_counters = [EventCounterCls(event_type, self.log_info) for
      event_type in EVENT_TYPES]
    #self.output_file = None
  #} end def

  def __del__(self): #{
    #if (None != self.output_file): #{
    #  self.output_file.Close()
    #} end if
    CloseLogFile(self)
  #} end def

  def SummarizePredictions(self): #{
    LogMsg(self, "Summarizing Barnacle predictions...")
    start = time.time()
    # initialize the library iterator for loading predictions
    lib_iterator = LibIteratorCls(self.options.input_path,
       self.LoadLibrary, self.options, self.log_info)
    # do not try to find event paths
    lib_iterator.options.get_paths = False
    # load the predictions from each library
    lib_iterator.IterateOverAllLibs()
    # create output file (OR JUST USE LOG FILE??)
    #self.CreateOutputFile()
    for event_counter in self.event_counters: #{
      #DebugMsg(self, "%s:\n%s" % (event_counter.key,
      #  "\n".join(["%s: %s -- %s" % (name, event_counter.gene_keys[name],
      #  ",".join(sorted(event_counter.genes[event_counter.gene_keys[name]])))
      #  for name in sorted(event_counter.gene_keys)])))
      # amalgamate gene keys
      event_counter.AmalgamateGeneKeys()
      #DebugMsg(self, "%s:\n%s" % (event_counter.key,
      #  "\n".join(["%s: %s -- %s" % (name, event_counter.gene_keys[name],
      #  ",".join(sorted(event_counter.genes[event_counter.gene_keys[name]])))
      #  for name in sorted(event_counter.gene_keys)])))
      # summarize event type
      #raise PredictionSummarizerError("not implemented!")
      self.SummarizeEventType(event_counter)
    #} end for
    if (self.options.write_events): #{
      # initialize the library iterator for writing recurrent predictions
      lib_iterator = LibIteratorCls(self.options.input_path,
         self.WriteEvents, self.options, self.log_info)
      # do not try to find event paths
      lib_iterator.options.get_paths = False
      # write the recurrent predictions for each library
      lib_iterator.IterateOverAllLibs()
      DebugMsg(self, "Num libs processed: %i" % lib_iterator.num_libs)
    #} end if
    LogMsg(self, "Total time spent summarizing Barnacle predictions: %s" %
      TimeSpent(start))
  #} end def

  def LoadLibrary(self, lib_info): #{
    start = time.time()
    LogMsg(self, "Loading predictions for library %s, output version: %s..." %
      (lib_info.lib_name, lib_info.barnacle_ver))
    for event_counter in self.event_counters: #{
      #DebugMsg(self, "Loading %s predictions..." % event_counter.key)
      event_counter.AddLibrary(lib_info.lib_name)
      event_file_name = "%s.barnacle.%s" % (lib_info.lib_name,
        event_counter.ext)
      event_path = os.path.join(lib_info.lib_dir, "8_predicted_events",
        event_file_name)
      #DebugMsg(self, "PATH: %s" % event_path)
      if (not os.path.isfile(event_path)): #{
        DebugMsg(self, "No %s predictions for library." % event_counter.key)
        continue
      #} end if
      event_parser = CandidateGroupParserCls(event_path)
      for event in event_parser: #{
        #DebugMsg(self, "Loading event %i" % event.id)
        event_counter.LoadEvent(lib_info.lib_name, event)
      #} end for
    #} end for
    LogMsg(self, "Time loading predictions: %s" % TimeSpent(start))
  #} end def

  #def CreateOutputFile(self): #{
  #  input_file_name = os.path.basename(self.options.input_path)
  #  input_root = os.path.splitext(input_file_name)[0]
  #  output_file_name = "%s.barnacle_summary.txt" % input_root
  #  output_path = os.path.join(self.options.libs_dir, output_file_name)
  #  self.output_file = FileBoxCls(output_path, "w",
  #    "cannot create output path")
  #} end def

  def SummarizeEventType(self, event_counter): #{
    start = time.time()
    LogMsg(self, "Summarizing %s predictions..." % event_counter.key)
    total_gene_set = set()
    total_gene_set_rec = set()
    for lib_key in sorted(event_counter.libraries.keys()): #{
      num_predictions = len(event_counter.libraries[lib_key])
      num_recurrent_predictions = 0
      gene_set = set()
      recurrent_gene_set = set()
      for gene_names_str in event_counter.libraries[lib_key]: #{
        gene_key = event_counter.gene_keys[gene_names_str]
        gene_set.add(gene_key)
        total_gene_set.add(gene_key)
        if (1 < len(event_counter.genes[gene_key])): #{
          num_recurrent_predictions += 1
          recurrent_gene_set.add(gene_key)
          total_gene_set_rec.add(gene_key)
        #} end if
      #} end for
      num_genes = len(gene_set)
      num_recurrent_genes = len(recurrent_gene_set)
      LogMsg(self, "%s: %i predictions (%i recurrent); %i genes "
        "(%i recurrent)" % (lib_key, num_predictions,
        num_recurrent_predictions, num_genes, num_recurrent_genes))
    #} end for
    LogMsg(self, "Total genes: %i (%i recurrent)" % (len(total_gene_set),
      len(total_gene_set_rec)))
    LogMsg(self, "Time summarizing predictions: %s" % TimeSpent(start))
  #} end def

  def WriteEvents(self, lib_info): #{
    start = time.time()
    LogMsg(self, "Writing recurrent and library-specific predictions for "
      "library %s, output version: %s..." % (lib_info.lib_name,
      lib_info.barnacle_ver))
    for event_counter in self.event_counters: #{
      DebugMsg(self, "Writing %s predictions..." %
        event_counter.key)
      event_file_name = "%s.barnacle.%s" % (lib_info.lib_name,
        event_counter.ext)
      event_path = os.path.join(lib_info.lib_dir, "8_predicted_events",
        event_file_name)
      #DebugMsg(self, "PATH: %s" % event_path)
      if (not os.path.isfile(event_path)): #{
        DebugMsg(self, "No %s predictions for library." % event_counter.key)
        continue
      #} end if
      event_parser = CandidateGroupParserCls(event_path)
      recurrent_path = "%s.rec" % event_path
      recurrent_file = FileBoxCls(recurrent_path, "w",
        "cannot create recurrent prediction file")
      lib_specific_path = "%s.lsp" % event_path
      lib_specific_file = FileBoxCls(lib_specific_path, "w",
        "cannot create library-specific prediction file")
      for event in event_parser: #{
        # get the gene_names_str for the event
        gene_names_str = event_counter.GetGeneNamesString(event)
        # get the appropriate gene_key
        gene_key = event_counter.gene_keys[gene_names_str]
        # check the number of libraries with predictions
        num_libs = len(event_counter.genes[gene_key])
        DebugMsg(self, "G%i: %s = %s (%s) = %i" % (event.id,
          gene_names_str, gene_key, ",".join(event_counter.genes[gene_key]),
          num_libs))
        # write the event if it is recurrent
        if (1 < num_libs): #{
          recurrent_file.Write(event.FullDataString())
        else:
          lib_specific_file.Write(event.FullDataString())
        #} end if
      #} end for
    #} end for
    LogMsg(self, "Time loading predictions: %s" % TimeSpent(start))
  #} end def
#} end class

class EventCounterCls: #{
  def __init__(self, event_info, log_info): #{
    self.key = event_info[0]
    self.ext = event_info[1]
    # libraries[lib_name] = list of gene_names strings
    self.libraries = dict()
    # gene_keys[gene_names string] = gene_key
    self.gene_keys = dict()
    # genes[gene_key] = lib_set
    self.genes = dict()
    self.log_info = log_info
  #} end def

  def AddLibrary(self, lib_key): #{
    self.libraries[lib_key] = list()
  #} end def

  def LoadEvent(self, lib_key, event): #{
    gene_names_str = self.GetGeneNamesString(event)
    self.libraries[lib_key].append(gene_names_str)
    gene_key = gene_names_str.lower()
    #DebugMsg(self, "Genes: %s (key: %s)" % (gene_names_str, gene_key))
    if (gene_names_str in self.gene_keys and
        gene_key != self.gene_keys[gene_names_str]): #{
      raise PredictionSummarizerError("inconsistent keys for %s: %s and %s" %
        (gene_names_str, gene_key, self.gene_keys[gene_names_str]))
    #} end if
    self.gene_keys[gene_names_str] = gene_key
    if (gene_key not in self.genes): #{
      self.genes[gene_key] = set()
    #} end if
    self.genes[gene_key].add(lib_key)
  #} end def

  def GetGeneNamesString(self, event): #{
    if ("fusion" in self.key): #{
      gene_sets = self.GetGeneSets(event)
      gene_sets[0].discard("no_exons")
      gene_sets[1].discard("no_exons")
      gene_names_str_1 = ",".join(sorted(gene_sets[0]))
      gene_names_str_2 = ",".join(sorted(gene_sets[1]))
      return ";".join(sorted([gene_names_str_1, gene_names_str_2]))
    else:
      gene_set = self.GetGeneSet(event)
      gene_set.discard("no_exons")
      return ",".join(sorted(gene_set))
    #} end if
  #} end def

  def GetGeneSets(self, event): #{
    gene_set_1 = set(event.members[0].genes_A)
    gene_set_2 = set(event.members[0].genes_B)
    #DebugMsg(self, "Initializing gene sets with: %s" %
    #  event.members[0].MainGenesString())
    if (1 == len(event.members)): #{
      return (gene_set_1, gene_set_2)
    #} end if
    for member in event.members[1:]: #{
      #DebugMsg(self, "Updating gene sets with: %s" %
      #  member.MainGenesString())
      overlap = (len(gene_set_1.intersection(member.genes_A)) +
        len(gene_set_2.intersection(member.genes_B)))
      swap_overlap = (len(gene_set_1.intersection(member.genes_B)) +
        len(gene_set_2.intersection(member.genes_A)))
      #DebugMsg(self, "Overlap: %i, Swap overlap: %i" % (overlap, swap_overlap))
      if (overlap > swap_overlap): #{
        #DebugMsg(self, "Not swapping")
        gene_set_1.update(member.genes_A)
        gene_set_2.update(member.genes_B)
      else:
        #DebugMsg(self, "Swapping")
        gene_set_1.update(member.genes_B)
        gene_set_2.update(member.genes_A)
      #} end if
    #} end for
    return (gene_set_1, gene_set_2)
  #} end def

  def GetGeneSet(self, event): #{
    gene_set = set()
    for member in event.members: #{
      gene_set.update(member.genes_A)
    #} end for
    return gene_set
  #} end def

  def AmalgamateGeneKeys(self): #{
    start = time.time()
    LogMsg(self, "Amalgamating %s gene keys..." % self.key)
    processed = set()
    for gene_names_str in sorted(self.gene_keys, key=lambda x: x.count(","),
        reverse=True): #{
      full_key = self.gene_keys[gene_names_str]
      DebugMsg(self, "%s: %s" % (gene_names_str, full_key))
      if (gene_names_str in processed): #{
        DebugMsg(self, "Skipping processed key")
        continue
      #} end if
      if (0 == gene_names_str.count(",")): #{
        DebugMsg(self, "Done")
        break
      #} end if
      if ("fusion" in self.key): #{
        self.AmalgamateGenePairKeys(gene_names_str, full_key, processed)
      else:
        for subset in PowerSet(gene_names_str.split(",")): #{
          DebugMsg(self, "Subset: %s" % ",".join(sorted(subset)))
          if (0 == len(subset) or
              len(gene_names_str.split(",")) == len(subset)): #{
            continue
          #} end if
          self.CheckSubKey(subset, None, full_key, processed)
        #} end for
      #} end if
    #} end for
    LogMsg(self, "Time amalgamating gene keys: %s" % TimeSpent(start))
  #} end def

  def AmalgamateGenePairKeys(self, gene_names_str, full_key, processed): #{
    (genes_1, genes_2) = sorted(gene_names_str.split(";"),
      key=lambda x: x.count(","), reverse=True)
    for subset_1 in PowerSet(genes_1.split(",")): #{
      #DebugMsg(self, "Subset1: %s" % ",".join(sorted(subset_1)))
      if (0 == len(subset_1)): #{
        continue
      #} end if
      #if ("," not in genes_2): #{
      #  self.CheckSubKey(subset_1, [genes_2], processed)
      #  continue
      #} end if
      for subset_2 in PowerSet(genes_2.split(",")): #{
        #DebugMsg(self, "Subset2: %s" % ",".join(sorted(subset_2)))
        if (0 == len(subset_2) or
            (len(genes_1.split(",")) == len(subset_1) and
             len(genes_2.split(",")) == len(subset_2))): #{
          continue
        #} end if
        self.CheckSubKey(subset_1, subset_2, full_key, processed)
      #} end for
    #} end for
  #} end def

  def CheckSubKey(self, subset_1, subset_2, full_key, processed): #{
    if ("fusion" in self.key): #{
      sub_names_str_1 = ",".join(sorted(subset_1))
      sub_names_str_2 = ",".join(sorted(subset_2))
      sub_names_str = ";".join(sorted([sub_names_str_1, sub_names_str_2]))
    else:
      sub_names_str = ",".join(sorted(subset_1))
    #} end if
    if (sub_names_str in processed): #{
      DebugMsg(self, "WARNING: sub key already processed: %s" % sub_names_str)
      return
    #} end if
    processed.add(sub_names_str)
    sub_key = sub_names_str.lower()
    DebugMsg(self, "SUBKEY %s: %s" % (sub_names_str, sub_key))
    if (sub_names_str not in self.gene_keys): #{
      DebugMsg(self, "Sub Key not present")
      return
    #} end if
    DebugMsg(self, "Amalgamating sub key with %s" % full_key)
    self.gene_keys[sub_names_str] = full_key
    if (sub_key in self.genes): #{
      self.genes[full_key].update(self.genes[sub_key])
      del self.genes[sub_key]
    #} end if
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class PredictionSummarizerError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Counts number of predictions, number of recurrent "
    "predictions, number of genes with predictions, and number of genes with "
    "recurrent predictions for each library in the input list.")
  args = [ "LIBRARY_LIST", "LIBS_DIR", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--write-events",
                    action="store_true", dest="write_events",
                    help="Write events to recurrent or library-specific "
                      "output files. [default]")
  parser.add_option("--count-only",
                    action="store_false", dest="write_events",
                    help="Do not write events to recurrent or "
                      "library-specific output files.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.set_defaults(write_events=True,
                      force=False,
                      debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (not opts_good): #{
    ErrMsg("bad option") #TODO
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.input_path, "libraries list", path_errors)
  # check the libraries path
  CheckDirPath(options.libs_dir, "libraries", path_errors,
    create=False)
  # get the log file name
  options.log_file_name = GetLogPath(options.input_path,
    "summary", options.libs_dir)
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
  #} end if
  # the paths are good if there are no path errors and no conflicting options
  return (opts_good and 0 == len(path_errors))
#} end def

def Main(): #{
  # get options and arguments
  parser = SetupOptionsParser()
  (options, args) = parser.parse_args()
  # if the right number of args was used
  if (parser.num_args == len(args)): #{
    options.input_path = EnsureAbsPath(args[0])
    options.libs_dir   = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      try: #{
        main_class_object = PredictionSummarizerCls(options)
        WriteCommand(main_class_object, sys.argv)
        main_class_object.SummarizePredictions()
      except (MyError), e:
        ErrMsg("ERROR while summarizing Barnacle predictions:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify the path to a file containing a list of "
      "libraries with Barnacle results (LIBRARY_LIST); and the path to the "
      "base directory containing the data for those libraries (LIBS_DIR).")
    return ES_OPT_ERR
  #} end if
  return ES_SUCCESS
#} end def

if __name__ == '__main__': #{
  try: #{
    exit_status = Main()
  except Exception, e:
    traceback.print_exc()
    exit_status = ES_EXCEPTION
  #} end try
  if (ES_SUCCESS == exit_status): #{
    print MSG_SUCCESS
  else:
    print MSG_FAIL
  #} end if
  sys.exit(exit_status)
#} end if
