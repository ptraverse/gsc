#! /usr/bin/env python
"""
get_recurrence.py

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
from utils.general import (SetupMainClass, TimeSpent, WriteCommand)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import CheckFilePath, EnsureAbsPath

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"

class RecurrenceCounterCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    # genes_in_lib[event_type][lib_name][gene_id] = gene_name
    self.genes_in_lib = dict([(predictor.key, dict()) for
      predictor in self.predictors])
    self.genes_in_lib['fusion_pairs'] = dict()
    # genes_with_events[event_type][gene_id] = gene_name
    self.genes_with_events = dict([(predictor.key, dict()) for
      predictor in self.predictors])
    self.genes_with_events['fusion_pairs'] = dict()
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def MainFunction(self): #{
    LogMsg(self, "RUNNING MAIN FUNCTION...") #TODO fix message
    start = time.time()
    #TODO implement the main function
    if (self.options.list_input): #{
      lib_iterator = LibIteratorCls(self.options.barnacle_path,
         self.GetBioTypesInFile, self.options, self.log_info)
      lib_iterator.list_of_paths = True
      lib_iterator.IterateOverAllLibs()
      DebugMsg(self, "Num libs found: %i" % lib_iterator.num_libs)
    else:
      lib_info = LibraryInfoCls(self.options, self.log_file)
      input_file_name = os.path.basename(self.options.barnacle_path)
      input_root = os.path.splitext(input_file_name)[0]
      lib_info.lib_name = input_root.replace(".anomalous_contig", "")
      lib_info.event_paths = [self.options.barnacle_path]
      self.GetBioTypesInFile(lib_info)
    #} end if
    LogMsg(self,
      "Time spent RUNNING MAIN FUNCTION: %s" % TimeSpent(start)) #TODO fix message
  #} end def

  def GetBioTypesInFile(self, lib_info): #{
    start = time.time()
    LogMsg(self, "Getting biologically typed events in file...")
    if (0 == len(lib_info.event_paths)): #{
      LogMsg(self, "Warning: no event paths found for library: %s\n" % (lib_info.lib_name))
      return
    #} end if
    if (1 < len(lib_info.event_paths)): #{
      LogMsg(self, "Warning: too many event paths found for library: %s\n  " % (lib_info.lib_name) + "\n  ".join(lib_info.event_paths))
    #} end if
    LogMsg(self, "  %s" % lib_info.event_paths[0])
    for predictor in self.predictors: #{
      self.genes_in_lib[predictor.key][lib_info.lib_name] = dict()
    #} end for
    self.genes_in_lib['fusion_pairs'][lib_info.lib_name] = dict()
    event_parser = CandidateGroupParserCls(lib_info.event_paths[0])
    output_files = self.CreateOutputFiles(lib_info.event_paths[0])
    for output_file in output_files.itervalues(): #{
      output_file.Close()
    #} end for
    self.OutputCounts()
    LogMsg(self, "Time processing file: %s" % TimeSpent(start))
  #} end def

  def CreateOutputFiles(self, barnacle_path): #{
    input_file_name = os.path.basename(barnacle_path)
    input_root = os.path.splitext(input_file_name)[0]
    output_files = dict()
    # setup the coordinates file for rechecking exon overlaps
    self.SetupEventCoordsFile(input_root, output_files)
    for predictor in self.predictors: #{
      output_file_name = "%s.%s" % (input_root, predictor.ext)
      output_path = os.path.join(self.options.output_dir, output_file_name)
      output_files[predictor.key] = FileBoxCls(output_path, "w",
        "cannot create %s output file" % predictor.description)
      if (self.options.pretty): #{
        pretty_path = "%s.pretty" % output_path
        output_files["%s_pretty" % predictor.key] = FileBoxCls(pretty_path,
          "w", "cannot create pretty %s output file" % predictor.description)
      #} end if
    #} end for
    return output_files
  #} end def

  def GetGenePairs(self, event): #{
    DebugMsg(self, "Getting gene pairs for event %s..." % event.id)
    gene_pairs = set()
    for member in event.members: #{
      if (0 == len(member.genes_A) or
          0 == len(member.genes_B)):
        continue
      #} end if
      for gene_A in member.genes_A: #{
        DebugMsg(self, "Gene A: %s" % gene_A)
        gene_A = RemoveExonTypeLabel(gene_A)
        if ("no_exons" == gene_A.lower()): #{
          continue
        #} end if
        for gene_B in member.genes_B: #{
          if ("no_exons" == gene_B.lower() or
              gene_A.lower() == gene_B.lower()):
            continue
          #} end if
          gene_B = RemoveExonTypeLabel(gene_B)
          gene_pair = [gene_A, gene_B]
          gene_pair_str = "/".join(sorted(gene_pair))
          gene_pairs.add(gene_pair_str)
        #} end for
      #} end for
    #} end for
    return gene_pairs
  #} end def

  def ProcessBiologicalEvent(self, type_key, lib_name, event, output_files): #{
    # update the counts
    new_genes = event.GeneSet(retain_case=True)
    if ("itds" == type_key): #{
      # remove UTR and non_coding genes
      new_genes = filter(IsNotUTRorNonCoding, new_genes)
    #} end if
    self.UpdateGeneList(self.genes_in_lib[type_key][lib_name], new_genes)
    self.UpdateGeneList(self.genes_with_events[type_key], new_genes)
    if ("fusions" == type_key): #{
      gene_pairs = self.GetGenePairs(event)
      self.UpdateGeneList(self.genes_in_lib['fusion_pairs'][lib_name],
        gene_pairs)
      self.UpdateGeneList(self.genes_with_events['fusion_pairs'], gene_pairs)
    #} end if
  #} end def

  def UpdateGeneList(self, gene_list, new_genes): #{
    for gene_name in new_genes: #{
      gene_id = gene_name.lower()
      if ("no_exons" != gene_id and
          gene_id not in gene_list):
        gene_list[gene_id] = gene_name
      #} end if
    #} end for
  #} end def

  def OutputCounts(self): #{
    for predictor in self.predictors: #{
      counts_path = os.path.join(self.options.output_dir,
        "%s.counts" % predictor.key)
      counts_file = FileBoxCls(counts_path, "w",
        "cannot create %s counts file" % predictor.description)
      self.WriteCount(counts_file, "Number of genes with %s events" %
        predictor.description.rstrip("s"))
      for lib in sorted(self.genes_in_lib[predictor.key].keys()): #{
        # remove "no_exons"
        #self.genes_in_lib[predictor.key][lib].discard("no_exons")
        self.WriteCount(counts_file, "%s: %i" %
          (lib, len(self.genes_in_lib[predictor.key][lib])))
        msg = "  Genes: "
        if (0 == len(self.genes_in_lib[predictor.key][lib])): #{
          msg += "None"
        else:
          msg += ", ".join(
            sorted(self.genes_in_lib[predictor.key][lib].itervalues()))
        #} end if
        self.WriteCount(counts_file, msg)
        if ("fusions" == predictor.key): #{
          msg = ("  Gene pairs (%i): " %
            len(self.genes_in_lib['fusion_pairs'][lib]))
          if (0 == len(self.genes_in_lib['fusion_pairs'][lib])): #{
            msg += "None"
          else:
            msg += ", ".join(
              sorted(self.genes_in_lib['fusion_pairs'][lib].itervalues()))
          #} end if
          self.WriteCount(counts_file, msg)
        #} end for
      #} end for
      self.WriteCount(counts_file, "Across all libraries: %i genes " %
        len(self.genes_with_events[predictor.key]) +
        "with at least one %s event in at least one library." %
        predictor.description.rstrip("s"))
      self.WriteCount(counts_file, "Genes:\n%s" %
        ", ".join(sorted(self.genes_with_events[predictor.key].itervalues())))
      if ("fusions" == predictor.key): #{
        self.WriteCount(counts_file, "Across all libraries: %i gene pairs " %
          len(self.genes_with_events['fusion_pairs']) +
          "with at least one %s event in at least one library." %
          predictor.description.rstrip("s"))
        self.WriteCount(counts_file, "Gene pairs:\n%s" % ", ".join(
          sorted(self.genes_with_events['fusion_pairs'].itervalues())))
      #} end if
      counts_file.Close()
    #} end for
  #} end def

  def WriteCount(self, counts_file, msg): #{
    counts_file.Write("%s\n" % msg)
    DebugMsg(self, msg)
  #} end def
#} end class

class EventAndGeneSetsCls: #{
  def __init__(self, event): #{
    self.event = event
    # a member-wise double-dictionary of gene sets
    # two gene sets per member (A and B)
    self.gene_set_dict = dict()
    for member in event.members: #{
      self.gene_set_dict[member.member_id] = {'A': set(), 'B': set()}
    #} end def
  #} end def

  def AddExons(self, member_id, region_id, new_exons): #{
    self.gene_set_dict[member_id][region_id].update(new_exons)
  #} end def

  def GetExons(self, member_id, region_id): #{
    return self.gene_set_dict[member_id][region_id]
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class RecurrenceCounterError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("") #TODO
  args = [ "LIB", "BARNACLE_FILE", ] #TODO
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  misc_group = OptionGroup(parser, "Miscellaneous Options")
  misc_group.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, even if the output "
                         "directory already exists.")
  misc_group.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  misc_group.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.add_option_group(misc_group)
  parser.set_defaults(force=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (not opts_good): #{
    ErrMsg("bad option") #TODO
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.barnacle_path, "Barnacle predictions", path_errors) #TODO
  # get and check the output path
  options.output_dir = GetOutDir(os.path.dirname(options.barnacle_path),
    "TASK_DESC") #TODO
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors,
      create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.barnacle_path,
      "TASK_DESC", options.output_dir) #TODO
  #} end if
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
    options.lib           = args[0]
    options.barnacle_path = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      try: #{
        main_class_object = RecurrenceCounterCls(options)
        WriteCommand(main_class_object, sys.argv)
        main_class_object.MainFunction()
      except (MyError), e:
        ErrMsg("ERROR while RUNNING MAIN FUNCTION:\n  %s" % e) #TODO
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify a library name (LIB); and the path to a "
      "Barnacle data file for that library (BARNACLE_FILE).")
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
