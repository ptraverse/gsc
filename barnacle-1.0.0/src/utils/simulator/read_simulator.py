#! /usr/bin/env python
"""
read_simulator.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, random, re, shlex, shutil, sys, time, traceback

# ensure that the sys.path is set up to allow importing Barnacle modules
barnacle_dir = sys.path[0]
if ("" == barnacle_dir): #{
  barnacle_dir = os.getcwd()
#} end if
while (not os.path.isfile(os.path.join(barnacle_dir, "barnacle.pl"))): #{
  barnacle_dir = os.path.dirname(barnacle_dir)
  if ("" == barnacle_dir or "/" == barnacle_dir): #{
    print "cannot find BARNACLE directory"
    raise Exception
  #} end if
#} end while
barnacle_dir = os.path.expanduser(barnacle_dir)
barnacle_dir = os.path.abspath(barnacle_dir)
if (barnacle_dir not in sys.path): #{
  sys.path.insert(0, barnacle_dir)
#} end if

# import custom modules
from version import VERSION
from utils.log import GetLogPath, CloseLogFile
from utils.error import MyError
from utils.general import (SetupMainClass, TimeSpent, WriteCommand,
    CheckConfigCommands, GetCommand, NormalizeChrID, IntFloor, IntCeiling)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, CheckNewFilePath, EnsureDirectoryExists, FileBoxCls, CleanLine)
from utils.subprocesses import RunCommandFromList, STREAM_OUT#, STRING_IN, STRING_OUT
from parsers.fasta import FastaFileCls
from parsers.fastq import FastqFileCls
from parsers.tokenizer import TokenizerCls, GetFieldAndValue
from parsers.genes.annotation import GeneAnnotationParserCls
from utils.multi_dict import SetupMultiDict, CheckMultiDict

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "SUCCESS"
MSG_FAIL = "FAIL"
RETIRED = "<retired>"
ID_UPDATE_HEADER = "Old stable ID, New stable ID, Release, Mapping score"
ST_MISSING = "missing"
ST_RETIRED = "retired"
ST_UPDATED = "updated"

class MainClassCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
    #CheckConfigCommands(self, ["dwgsim", "idmapper", "idmapper_plus"])
    #CheckConfigCommands(self, ["dwgsim"])
    self.use_gene_id = True
    # gnpairs[gene_id] = total number of read-pairs for all isoforms of gene
    self.gnpairs = dict()
    # tnpairs[gene_id][transcript_id] =
    #   number of read-pairs to generate for each isoform
    self.tnpairs = dict()
    # t2g_map[transcript_id] = gene_id
    self.t2g_map = dict()
    self.read_files = dict()
    self.gid_updates = dict()
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
    if (hasattr(self, "read_files")): #{
      for file in self.read_files.itervalues(): #{
        file.Close()
      #} end for
    #} end if
  #} end def

  def Run(self): #{
    LogMsg(self, "Simulating reads from wildtype and event sequences...")
    start = time.time()
    self.ConsolidateIsoforms()
    self.LoadWildTypeCoverages()
    self.OpenReadFiles()
    self.GenerateWildtypeReads()
    self.GenerateEventReads()
    LogMsg(self, "Time spent simulating reads from wildtype and event sequences: %s" %
      TimeSpent(start))
  #} end def

  def ConsolidateIsoforms(self): #{
    LogMsg(self, "Consolidating gene isoforms...")
    start = time.time()
    annot_file = GeneAnnotationParserCls(self.options.annot_path,
        log_info=self.log_info)
    for transcript in annot_file: #{
      if (None != self.options.chr_filter and
          NormalizeChrID(transcript.chrom) not in self.options.chr_filter): #{
        continue
      #} end if
      #ExtremeDebugMsg(self, "  transcript: %s" % transcript)
      tid = transcript.transcript_id
      if (hasattr(transcript, "gene_id")): #{
        gene_id = transcript.gene_id
      else:
        gene_id = transcript.gene_name
        self.use_gene_id = False
      #} end if
      SetupMultiDict(self.tnpairs, [gene_id, tid],
          0, list_end=False)
      if (tid in self.t2g_map and
          self.t2g_map[tid] != gene_id): #{
        raise MainClassError("Conflicting gene ids for transcript %s:"
            " %s != %s" % (tid, self.t2g_map[tid], gene_id))
      #} end if
      self.t2g_map[tid] = gene_id
    #} end for
    annot_file.Close()
    if (self.log_info["extreme_debug"]): #{
      for gene_id in sorted(self.tnpairs.keys()): #{
        LogMsg(self, "  Gene: %s" % gene_id)
        for tid in sorted(self.tnpairs[gene_id].keys()): #{
          LogMsg(self, "    Transcript: %s" % tid)
        #} end for
      #} end for
    #} end if
    LogMsg(self, "Time spent consolidating gene isoforms: %s" %
      TimeSpent(start))
  #} end def

  def LoadWildTypeCoverages(self): #{
    LogMsg(self, "Loading wildtype coverages...")
    start = time.time()
    cov_file = FileBoxCls(self.options.wcov_path, "r",
        "cannot read wildtype coverages file")
    num_no_reads = 0
    num_genes = 0
    num_missing = 0
    missing_path = os.path.join(self.options.woutput_dir, "missing.txt")
    missing_file = FileBoxCls(missing_path, "w",
        "cannot create missing wildtype genes file")
    missing_coverage = dict()
    for line in cov_file: #{
      gene_coverage = CoverageCls(line, log_info=self.log_info)
      if (None != self.options.chr_filter and
          gene_coverage.chrom not in self.options.chr_filter): #{
        continue
      #} end if
      msg = "Loading coverage for %s (%s)" % (gene_coverage.gene_name,
          gene_coverage.gene_id)
      ExtremeDebugMsg(self, "  %s" % msg)
      if (1 > gene_coverage.exon_nreads): #{
        ExtremeDebugMsg(self, "    Coverage record with no reads: %f" %
            gene_coverage.exon_nreads)
        num_no_reads += 1
        continue
      #} end if
      num_genes += 1
      self.tracking = list()
      self.tracking.append(msg)
      if (self.use_gene_id): #{
        gene_id = gene_coverage.gene_id
      else:
        gene_id = gene_coverage.gene_name
      #} end if
      if (gene_id not in self.tnpairs): #{
        ExtremeDebugMsg(self, "    Gene %s not in annotations" % gene_id)
        if (gene_id in missing_coverage): #{
          raise MainClassError("Duplicated gene %s in wildtype coverages "
              "file!" % gene_id)
        #} end if
        missing_coverage[gene_id] = gene_coverage
        missing_file.WriteLine(gene_id)
        num_missing += 1
        continue
      #} end if
      self.SpreadReadsAmongstIsoforms(gene_id, gene_coverage.exon_nreads)
    #} end for
    if (0 < num_no_reads): #{
      LogMsg(self, "Coverage file has %i records with less than one read" %
          num_no_reads)
    #} end if
    if (0 < num_missing): #{
      LogMsg(self, "%i/%i gene IDs are not present in annotations file")
      #LogMsg(self, "%i/%i gene IDs are not present in annotations file, "
      #    "will attempt to update them..." % (num_missing, num_genes))
      #self.UpdateGeneIDs(missing_coverage, missing_file, plus=False)
    #} end if
    missing_file.Close()
    cov_file.Close()
    LogMsg(self, "Time spent loading wildtype coverages: %s" %
      TimeSpent(start))
  #} end def

  def SpreadReadsAmongstIsoforms(self, gene_id, num_reads): #{
    pairs_remaining = IntFloor(float(num_reads)/2.0)
    self.tracking.append("Total number of read-pairs for gene: %i" %
        pairs_remaining)
    if (gene_id not in self.gnpairs): #{
      self.gnpairs[gene_id] = 0
    #} end if
    self.gnpairs[gene_id] += pairs_remaining
    iso_keys = self.tnpairs[gene_id].keys()
    num_isoforms = len(iso_keys)
    if (1 < num_isoforms): #{
      #last_index = random.randint(0,num_isoforms-1)
      #isoforms = (iso_keys[:last_index] + iso_keys[last_index+1:])
      #random.shuffle(isoforms)
      random.shuffle(iso_keys)
    #else:
    #  last_index = 0
    #  isoforms = []
    #} end if
    num_expressed = 0
    #for tid in isoforms: #{}
    #for tid in iso_keys[:-1]: #{}
    for tid in iso_keys: #{
      # if isoform_fract > 1.0 - 1/n, where n is the number of isoforms
      # remaining then just give all the reads to the current isoform
      isoform_fract = random.random()
      threshold = 1.0 - (1.0/float(num_isoforms-num_expressed))
      msg = "Initial fract: %.3f, Threshold: %.3f" % (isoform_fract,
          threshold)
      ExtremeDebugMsg(self, "    %s" % msg)
      self.tracking.append(msg)
      if (isoform_fract > threshold): #{
        isoform_fract = 1.0
      #} end if
      isoform_npairs = IntCeiling(pairs_remaining*isoform_fract)
      self.tnpairs[gene_id][tid] += isoform_npairs
      msg = "Assigning %.2f (%i) of remaining (%i) to %s" % (isoform_fract,
          isoform_npairs, pairs_remaining, tid)
      ExtremeDebugMsg(self, "    %s" % msg)
      self.tracking.append(msg)
      pairs_remaining -= isoform_npairs
      msg = "new remaining: %i" % pairs_remaining
      ExtremeDebugMsg(self, "      %s" % msg)
      self.tracking.append(msg)
      num_expressed += 1
      if (0 >= pairs_remaining): #{
        #num_expressed -= 1
        break
      #} end if
    #} end for
    #last_key = iso_keys[last_index]
    #last_key = iso_keys[-1]
    #ExtremeDebugMsg(self, "    Assigning all remaining (%i) reads to %s" %
    #    (pairs_remaining, last_key))
    #self.tnpairs[gene_id][last_key] = pairs_remaining
    msg = "Added reads to %i of %i isoforms" % (num_expressed, num_isoforms)
    ExtremeDebugMsg(self, "    %s" % msg)
    self.tracking.append(msg)
    sum_tnpairs = sum(self.tnpairs[gene_id].itervalues())
    if (sum_tnpairs != self.gnpairs[gene_id]): #{
      raise MainClassError("error while calculating isoform coverage for "
        "gene %s, sum of isoform read-pair counts: %i != total gene "
        "read-pair count: %i\n%s\nUPDATED GENE IDS:\n%s" % (gene_id,
        sum_tnpairs, self.gnpairs[gene_id], "\n".join(self.tracking),
        "\n".join(("%s:%s" % (old_gid, ", ".join(self.gid_updates[old_gid]))
        for old_gid in sorted(self.gid_updates.keys())))))
    #} end if
  #} end def

  def UpdateGeneIDs(self, missing_coverage, missing_file, plus=False): #{
    DebugMsg(self, "Updating gene IDs (plus = %s)..." % plus)
    missing_file.Close()
    start = time.time()
    update_stream = GeneIDUpdateStreamCls(missing_file.path,
      self.options, plus=plus, log_info=self.log_info)
    #ExtremeDebugMsg(self, "Constructed update stream object!")
    #MAX_LOOPS = 20
    #count = 0
    # update[old_gid] = set(new_gids)
    for update in update_stream: #{
      ExtremeDebugMsg(self, "Processing update from stream.")
      for old_gid in update: #{
        ExtremeDebugMsg(self, "  Adding mapping from %s to %s" %
            (old_gid, ", ".join(update[old_gid])))
        if (old_gid not in self.gid_updates): #{
          self.gid_updates[old_gid] = set()
        #} end if
        self.gid_updates[old_gid].update(update[old_gid])
      #} end for
      #count += 1
      #if (MAX_LOOPS < count): #{
      #  raise MainClassError("Too many loops while processing updates from stream!")
      #} end if
    #} end for
    LogMsg(self, "Time spent running IDmapper script: %s" %
      TimeSpent(start))
    missing_file.Open()
    still_missing = dict()
    for old_gid in missing_coverage: #{
      ExtremeDebugMsg(self, "Checking whether %s could be updated..." %
          old_gid)
      if (not hasattr(missing_coverage[old_gid], "status")): #{
        missing_coverage[old_gid].status = ST_MISSING
      #} end if
      new_id_set = set()
      if (old_gid in self.gid_updates): #{
        new_id_set = self.gid_updates[old_gid]
        if (RETIRED in new_id_set and
            ST_UPDATED != missing_coverage[old_gid].status): #{
          missing_coverage[old_gid].status = ST_RETIRED
        #} end if
        DebugMsg(self, "  Found %i new IDs: %s" %
            (len(new_id_set), ", ".join(sorted(new_id_set))))
        new_id_set.intersection_update(self.tnpairs.keys())
        DebugMsg(self, "    %i of them in the annotations file: %s" %
            (len(new_id_set), ", ".join(sorted(new_id_set))))
      #} end if
      if (0 < len(new_id_set)): #{
        missing_coverage[old_gid].status = ST_UPDATED
        new_id = random.choice(list(new_id_set))
        DebugMsg(self, "  Randomly selected new ID: %s" % new_id)
        self.tracking.append("Updated gene ID %s to %s" % (old_gid, new_id))
        self.SpreadReadsAmongstIsoforms(new_id,
            missing_coverage[old_gid].exon_nreads)
      else:
        DebugMsg(self, "  Could not update gene ID %s!" % old_gid)
        still_missing[old_gid] = missing_coverage[old_gid]
        missing_file.WriteLine(old_gid)
      #} end if
    #} end for
    missing_file.Close()
    LogMsg(self, "Updated gene IDs for %i coverage values." %
        (len(missing_coverage)-len(still_missing)))
    if (0 < len(still_missing)): #{
      if (plus): #{
        LogMsg(self, "Skipped %i coverage values for genes missing from "
            "annotations file (or retired)." % len(still_missing))
      else:
        LogMsg(self, "%i/%i gene IDs could not be updated directly, will "
            "attempt to update them via transcript ids..." %
            (len(still_missing), len(missing_coverage)))
        self.UpdateGeneIDs(still_missing, missing_file, plus=True)
      #} end if
    #} end if
    if (not plus): #{
      missing_file.Open()
      for old_gid in sorted(missing_coverage.keys()): #{
        missing_file.WriteLine(missing_coverage[old_gid])
      #} end for
      missing_file.Close()
    #} end if
  #} end def

  #def OldUpdateGeneIDs(self, missing_coverage, missing_file): #{
  #  #num_updated = 0
  #  #num_missing = 0
  #  #old_gid = gene_id
  #  missing_file.Close()
  #  start = time.time()
  #  update_stream = GeneIDUpdateStreamCls(missing_file.path,
  #    self.options, log_info=self.log_info)
  #    #self.options.cfg, log_info=self.log_info, dpt=self.options.dpt)
  #  #ExtremeDebugMsg(self, "Constructed update stream object!")
  #  #MAX_LOOPS = 20
  #  #count = 0
  #  for update in update_stream: #{
  #    ExtremeDebugMsg(self, "Processing update from stream.")
  #    for old_gid in update: #{
  #      ExtremeDebugMsg(self, "  Adding mapping from %s to %s" %
  #          (old_gid, ", ".join(update[old_gid])))
  #      if (old_gid not in self.gid_updates): #{
  #        self.gid_updates[old_gid] = set()
  #      #} end if
  #      self.gid_updates[old_gid].update(update[old_gid])
  #    #} end for
  #    #count += 1
  #    #if (MAX_LOOPS < count): #{
  #    #  raise MainClassError("Too many loops while processing updates from stream!")
  #    #} end if
  #  #} end for
  #  LogMsg(self, "Time spent running IDmapper script: %s" %
  #    TimeSpent(start))
  #  #status = dict()
  #  missing_file.Open()
  #  still_missing = dict()
  #  for old_gid in missing_coverage: #{
  #    ExtremeDebugMsg(self, "Checking whether %s could be updated..." %
  #        old_gid)
  #    missing_coverage[old_gid].status = "missing"
  #    new_id_set = set()
  #    if (old_gid in self.gid_updates): #{
  #      new_id_set = self.gid_updates[old_gid]
  #      if (RETIRED in new_id_set): #{
  #        #status[old_gid] = "retired"
  #        missing_coverage[old_gid].status = "retired"
  #      #} end if
  #      DebugMsg(self, "  Found %i new IDs: %s" %
  #          (len(new_id_set), ", ".join(sorted(new_id_set))))
  #      new_id_set.intersection_update(self.tnpairs.keys())
  #      DebugMsg(self, "    %i of them in the annotations file: %s" %
  #          (len(new_id_set), ", ".join(sorted(new_id_set))))
  #    #} end if
  #    if (0 < len(new_id_set)): #{
  #      missing_coverage[old_gid].status = "updated"
  #      new_id = random.choice(list(new_id_set))
  #      DebugMsg(self, "  Randomly selected new ID: %s" % new_id)
  #      self.tracking.append("Updated gene ID %s to %s" % (old_gid, new_id))
  #      #num_updated += 1
  #      self.SpreadReadsAmongstIsoforms(new_id,
  #          missing_coverage[old_gid].exon_nreads)
  #    else:
  #      DebugMsg(self, "  Could not update gene ID %s!" % old_gid)
  #      still_missing[old_gid] = missing_coverage[old_gid]
  #      missing_file.WriteLine(old_gid)
  #      #missing_file.WriteLine("\t".join([line, status]))
  #      #num_missing += 1
  #      #if (old_gid not in status): #{
  #      #  status[old_gid] = "missing"
  #      #} end if
  #      #status = "missing"
  #      #if (retired): #{
  #      #  status = "retired"
  #      #} end if
  #      #return (None, status)
  #      #continue
  #    #} end if
  #    #return (new_id, "updated")
  #  #} end for
  #  missing_file.Close()
  #  #(gene_id, status) = self.UpdateGeneID(gene_id)
  #  #if (None == gene_id): #{
  #  #  num_missing += 1
  #  #  continue
  #  #else:
  #  #  num_updated += 1
  #  #} end if
  #  LogMsg(self, "Updated gene IDs for %i coverage values." %
  #      (len(missing_coverage)-len(still_missing)))
  #  if (0 < len(still_missing)): #{
  #    LogMsg(self, "%i/%i gene IDs could not be updated directly, will "
  #        "attempt to update them via transcript ids..." %
  #        (len(still_missing), len(missing_coverage)))
  #    self.UpdateGeneIDsPlus(still_missing, missing_file)
  #  #} end if
  #  missing_file.Open()
  #  for old_gid in sorted(missing_coverage.keys()): #{
  #    missing_file.WriteLine(missing_coverage[old_gid])
  #  #} end for
  #  missing_file.Close()
  #} end def

  def UpdateGeneIDsPlus(self, missing_coverage, missing_file): #{
    LogMsg(self, "Updated gene IDs for %i coverage values." %
        (len(missing_coverage)-len(still_missing)))
    if (0 < len(still_missing)): #{
      LogMsg(self, "Skipped %i coverage values for genes missing from "
          "annotations file (or retired)." % len(still_missing))
    #} end if
  #} end def

  #def UpdateGeneID(self, old_id): #{
  #  command = [GetCommand(self, "idmapper"), "-s", "human"]
  #  DebugMsg(self, "Updating gene ID %s:\n  %s" %
  #      (old_id, " ".join(command)))
  #  #output = RunCommandFromList(command, stdin=STRING_IN, stdout=STRING_OUT,
  #  #    input=old_id, dpt=self.options.dpt)[0]
  #  DebugMsg(self, "OUTPUT:\n%s" % output)
  #  retired = False
  #  new_id_set = set()
  #  # skip the header
  #  for update in output.split("\n")[1:]: #{
  #    if ("" == update): #{
  #      continue
  #    #} end if
  #    DebugMsg(self, "  Processing update: %s" % update)
  #    new_id = update.split(", ")[1]
  #    if ("." in new_id): #{
  #      new_id = new_id.split(".")[0]
  #    #} end if
  #    DebugMsg(self, "    New ID: %s" % new_id)
  #    if ("<retired>" == new_id): #{
  #      retired = True
  #    elif (new_id != old_id): #{
  #      new_id_set.add(new_id)
  #    #} end if
  #  #} end for
  #  DebugMsg(self, "Found %i new IDs: %s" %
  #      (len(new_id_set), ", ".join(sorted(new_id_set))))
  #  new_id_set.intersection_update(self.tnpairs.keys())
  #  DebugMsg(self, "  %i of them in the annotations file: %s" %
  #      (len(new_id_set), ", ".join(sorted(new_id_set))))
  #  if (0 == len(new_id_set)): #{
  #    DebugMsg(self, "Could not update gene ID!")
  #    status = "missing"
  #    if (retired): #{
  #      status = "retired"
  #    #} end if
  #    return (None, status)
  #  #} end if
  #  new_id = random.choice(list(new_id_set))
  #  DebugMsg(self, "Randomly selected new ID: %s" % new_id)
  #  return (new_id, "updated")
  #} end def

  def OpenReadFiles(self): #{
    desc = {
        "wread1":"wildtype read 1",
        "wread2":"wildtype read 2",
        "eread1":"event read 1",
        "eread2":"event read 2",
    }
    for type in desc.keys(): #{
      self.read_files[type] = FileBoxCls(getattr(self.options,
        "%s_seqs_path" % type), "w", "cannot create %s sequences file" %
        desc[type])
      self.read_files[type].Open()
    #} end for
  #} end def

  def GenerateWildtypeReads(self): #{
    LogMsg(self, "Generating wildtype reads...")
    start = time.time()
    seq_file = FastaFileCls(self.options.wseq_path,
        "cannot read sequences file")
    npairs_file = FileBoxCls(self.options.wnreads_path, "w",
        "cannot create wildtype read-pair counts file")
    # number of sequences from which read-pairs were actually simulated
    nseqs_sim = 0
    for seq_obj in seq_file: #{
      npairs = self.GetWTTranscriptReadCount(seq_obj)
      if (0 == npairs): #{
        continue
      #} end if
      nseqs_sim += 1
      #cov_line = cov_file.next()
      #coverage = float(cov_line) + self.options.cov_adjust
      #num_reads = (len(seq_obj) / self.options.read_length) * coverage
      coverage = (npairs*2 *
          float(self.options.read_length) / float(len(seq_obj)))
      npairs_file.WriteLine("%s %i %f" % (seq_obj.id, npairs, coverage))
      self.SimulateReads(seq_obj, npairs, "w")
    #} end for
    npairs_file.Close()
    seq_file.Close()
    LogMsg(self, "Simulated reads from %i wildtype sequences" % nseqs_sim)
    LogMsg(self, "Time spent generating wildtype reads: %s" %
      TimeSpent(start))
  #} end def

  def GetWTTranscriptReadCount(self, seq_obj): #{
    if (len(seq_obj) <= self.options.frag_length): #{
      ExtremeDebugMsg(self, "Sequence %s shorter than fragment length: "
          "%i < %i" % (seq_obj.id, len(seq_obj), self.options.frag_length))
      return 0
    #} end if
    if (seq_obj.id not in self.t2g_map): #{
      ExtremeDebugMsg(self, "No gene information for transcript %s "
          "has been loaded" % seq_obj.id)
      return 0
    #} end if
    gene_id = self.t2g_map[seq_obj.id]
    if (0 == self.gnpairs.get(gene_id, 0)): #{
      ExtremeDebugMsg(self, "No read-pairs to generate for gene %s (%s)" %
          (gene_id, seq_obj.id))
      return 0
    #} end if
    if (seq_obj.id not in self.tnpairs[gene_id]): #{
      ExtremeDebugMsg(self, "No transcript information for transcript "
          "%s has been loaded" % seq_obj.id)
      return 0
    #} end if
    npairs = self.tnpairs[gene_id][seq_obj.id]
    if (0 == npairs): #{
      ExtremeDebugMsg(self, "No reads to generate for transcript %s" %
          seq_obj.id)
    #} end if
    return npairs
  #} end def

  def GenerateEventReads(self): #{
    LogMsg(self, "Generating event reads...")
    start = time.time()
    seq_file = FastaFileCls(self.options.eseq_path,
        "cannot read sequences file")
    npairs_file = FileBoxCls(self.options.enreads_path, "w",
        "cannot create event read counts file")
    cov_file = FileBoxCls(self.options.ecov_path, "r",
        "cannot read event coverages file")
    # number of sequences from which reads were actually simulated
    nseqs_sim = 0
    for seq_obj in seq_file: #{
      if (len(seq_obj) <= self.options.frag_length): #{
        LogMsg(self, "Sequence %s shorter than fragment length: "
            "%i < %i" % (seq_obj.id, len(seq_obj), self.options.frag_length))
        continue
      #} end if
      nseqs_sim += 1
      seq_obj.covered = False
      while (not seq_obj.covered): #{
        cov_line = cov_file.next()
        #coverage = float(cov_line) + self.options.cov_adjust
        coverage = float(cov_line)
        nreads = (coverage *
            (float(len(seq_obj)) / float(self.options.read_length)))
        npairs = IntFloor(float(nreads) / 2.0)
        if (1 > npairs): #{
          ExtremeDebugMsg(self, "    coverage %.3f too low, no reads." %
              coverage)
          continue
        #} end if
        #coverage = nreads * self.options.read_length / len(seq_obj)
        self.SimulateReads(seq_obj, npairs, "e")
      #} end while
      npairs_file.WriteLine("%s %i %f" % (seq_obj.id, npairs, coverage))
    #} end for
    cov_file.Close()
    npairs_file.Close()
    seq_file.Close()
    LogMsg(self, "Simulated reads from %i event sequences" % nseqs_sim)
    LogMsg(self, "Time spent generating event reads: %s" %
      TimeSpent(start))
  #} end def

  # type should be "w" for wildtype, or "e" for event
  def SimulateReads(self, seq_obj, num_pairs, type): #{
    tmp_dir = os.path.join(getattr(self.options, "%soutput_dir" % type),
        "tmp")
    EnsureDirectoryExists(tmp_dir)
    dwgsim_log_path = os.path.join(tmp_dir, "dwgsim.log")
    # write the source sequence to a temporary file
    tmp_in_path = os.path.join(tmp_dir, "tmp_ref.fa")
    tmp_in_file = FileBoxCls(tmp_in_path, "w",
        "cannot create temporary sequence file")
    tmp_in_file.WriteLine(seq_obj.Output())
    tmp_in_file.Close()
    # write the reads to a temporary file, before appending
    tmp_out_path = os.path.join(tmp_dir, "out")
    # construct the dwgsim command
    dwgsim_path = os.path.join(os.environ["BARNACLE_PATH"],
        "utils", "dwgsim")
    #command = [GetCommand(self, "dwgsim"),
    command = [dwgsim_path,
      "-d%i" % self.options.frag_length,
      "-N%i" % num_pairs,
      "-1%i" % self.options.read_length,
      "-2%i" % self.options.read_length,
      "-r%f" % self.options.mut_rate,
      "-y0",  # probability of a random DNA read
              # (randomly choose base for each position)
    ]
    if (None != self.options.err_rate): #{
      command.append("-e%f" % self.options.err_rate)
      command.append("-E%f" % self.options.err_rate)
    #} end if
    if (None != self.options.std_dev): #{
      command.append("-s%f" % self.options.std_dev)
    #} end if
    if (None != self.options.extra_dwgsim): #{
      command.extend(shlex.split(self.options.extra_dwgsim))
    #} end if
    command.extend([tmp_in_path, tmp_out_path])
    ExtremeDebugMsg(self, "  Running dwgsim: %s" % " ".join(command))
    es = RunCommandFromList(command, dpt=self.options.dpt,
        stderr=dwgsim_log_path)
    if (self.log_info["extreme_debug"]): #{
      dwgsim_log = FileBoxCls(dwgsim_log_path, "r")
      for line in dwgsim_log: #{
        LogMsg(self, "    %s" % line)
      #} end for
    #} end if
    if (0 > es): #{
      raise MainClassError("dwgsim command terminated by signal %i\n  %s" %
          (-es, " ".join(command)))
    elif (0 < es):
      raise MainClassError("dwgsim command failed with exit status: %i\n  %s" %
          (es, " ".join(command)))
    #} end if
    if ("e" == type): #{
      if (self.IsEventCovered(seq_obj, tmp_out_path)): #{
        seq_obj.covered = True
      else:
        return
      #} end if
    #} end if
    self.AppendNewReads(tmp_out_path, type, 1)
    self.AppendNewReads(tmp_out_path, type, 2)
  #} end def

  def IsEventCovered(self, event, tmp_prefix): #{
    tmp_reads_path = "%s.bwa.read1.fastq" % (tmp_prefix)
    tmp_file = FastqFileCls(tmp_reads_path,
        "cannot read temporary reads file")
    for read in tmp_file: #{
      id_parts = read.id.split("_")
      starts = map(int, id_parts[1:3])
      for i in range(2): #{
        if (self.ReadCoversEvent(starts[i], event)): #{
          return True
        #} end if
      #} end for
    #} end for
  #} end def

  def ReadCoversEvent(self, rstart, event): #{
    split_pos = None
    extra_len = 0
    tokenizer = TokenizerCls(event.extra, " ")
    for token in tokenizer: #{
      try: #{
        (field, value) = GetFieldAndValue(token)
      except ValueError:
        raise MainClassError("cannot parse read field: %s" % token)
      #} end try
      if ("SPLIT" == field): #{
        split_pos = int(value)
      elif ("EXTRA" == field):
        extra_len = len(value)
      #} end if
    #} end for
    if (None == split_pos): #{
      raise MainClassError("Could not get event sequence split position "
          "from: %s %s" % (event.id, event.extra))
    #} end if
    rend = rstart + self.options.read_length
    if (rstart < (split_pos - self.options.min_overlap) and
        rend   > (split_pos + self.options.min_overlap)): #{
      return True
    #} end if
    return False
  #} end def

  # type should be "w" for wildtype, or "e" for event
  def AppendNewReads(self, tmp_prefix, type, index): #{
    tmp_reads_path = "%s.bwa.read%i.fastq" % (tmp_prefix, index)
    tmp_file = FileBoxCls(tmp_reads_path, "r",
        "cannot read temporary reads file")
    tmp_file.Open()
    key = "%sread%i" % (type, index)
    shutil.copyfileobj(tmp_file.file, self.read_files[key].file)
    tmp_file.Close()
  #} end def
#} end class

field_types = {
    "gene_id":"str",          #
    "transcript_id":"str",    # or "merged_<gene_id>" if the analysis was performed in collapse mode
    "chrom":"str",            #
    "tstart":"int",           # relative to the positive strand, so start<end
    "tend":"int",             # relative to the positive strand, so start<end
    "strand":"str",           # "+" or "-"
    "exon_length":"int",      # sum of the length of all exons in this transcript or collapsed gene
    "intron_length":"int",    # sum of the length of all introns in this transcript or collapsed gene
    "exon_nreads":"float",    # number of fractional reads inside this merged_gene or transcript's exons (sum of the fraction of each read inside all exons)
    "exon_tcov":"int",        # total coverage across all the exons in this merged_gene or transcript (sum of the coverage depth at each base in all exons)
    "intron_nreads":"float",  # number of fractional reads inside this merged_gene or transcript's introns (sum of the fraction of each read inside all introns)
    "intron_tcov":"int",      # total coverage across all the introns in this merged_gene or transcript (sum of the coverage depth at each base in all introns)
    "exon_acov":"float",      # average coverage over all exons -- sum of the coverage depth at each base in all exons divided by the sum of the exon lengths)
    "full_acov":"float",      # average coverage over all introns and exons -- sum of the coverage depth at each base between the merged_gene or transcript's start and end divided by the number of bases between the gene's start and end
    "rpkm":"float",           # normalized coverage (RPKM) -- (number of fractional reads in all exons in this merged gene or transcript x 1000000000)/(NORM_TOTAL x sum of the length of all exons in this merged gene or transcript)
    "gene_name":"str",        # or gene id if symbol is unavailable
    "bio_type":"str",         # or "-" if biotype unavailable
    "gene_desc":"str",        # or "-" if description is unavailable
}

cov_fields = [ "gene_id", "transcript_id", "chrom", "tstart", "tend",
    "strand", "exon_length", "intron_length", "exon_nreads", "exon_tcov",
    "intron_nreads", "intron_tcov", "exon_acov", "full_acov", "rpkm",
    "gene_name", "bio_type", "gene_desc",
]

class CoverageCls: #{
  def __init__(self, line, log_info=None): #{
    #ExtremeDebugMsg(log_info, "  Constructing coverage object. (%i)" %
    #    len(cov_fields))
    self.line = line
    tokenizer = TokenizerCls(line)
    for field in cov_fields: #{
      #ExtremeDebugMsg(log_info, "    Loading value for field %s" % field)
      if ("str" == field_types[field]): #{
        setattr(self, field, tokenizer.next())
      elif ("int" == field_types[field]):
        setattr(self, field, int(tokenizer.next()))
      elif ("float" == field_types[field]):
        setattr(self, field, float(tokenizer.next()))
      else:
        raise CoverageParseError("unrecognized coverage field type: %s" %
            field_types[field])
      #} end if
    #} end for
    self.chrom = NormalizeChrID(self.chrom)
    if (self.tstart > self.tend): #{
      raise CoverageParseError("invalid transcript coordinates: %i > %i" %
          (self.tstart, self.tend))
    #} end if
    if (self.strand not in ["+", "-"]): #{
      raise CoverageParseError("invalid transcript strand: %s, must be "
          "\"+\" or \"-\"." % self.strand)
    #} end if
  #} end def

  def __str__(self): #{
    if (hasattr(self, "status")): #{
      return "\t".join([self.line, self.status])
    #} end if
    return self.line
  #} end def
#} end class

class GeneIDUpdateStreamCls: #{
  def __init__(self, old_id_path, options, plus=False, log_info=None): #{
    self.log_info = log_info
    self.options = options
    tool_name = "idmapper"
    if (plus): #{
      tool_name = "idmapper_plus"
    #} end if
    self.command = [GetCommand(self, tool_name), "-s", "human", "-f", old_id_path]
    DebugMsg(self, "Update command:\n  %s" % (" ".join(self.command)))
    self.stream = RunCommandFromList(self.command, stdout=STREAM_OUT,
        dpt=getattr(self.options, "dpt", False))
    #ExtremeDebugMsg(self, "  finished creating update stream")
    self.GetLine()
    #ExtremeDebugMsg(self, "  got the first line")
  #} end def

  def __iter__(self): #{
    #ExtremeDebugMsg(self, "Creating iterator for GeneIDUpdateStreamCls")
    return self
  #} end def

  # returns a dictionary: mapping[old_id] = set(new_ids)
  def next(self): #{
    #ExtremeDebugMsg(self, "In GeneIDUpdateStreamCls.next()...")
    if (getattr(self, "finished", False)): #{
      raise StopIteration
    #} end if
    if (ID_UPDATE_HEADER != self.curr_line): #{
      raise MainClassError("Invalid ID update header line: \"%s\"" %
          self.curr_line)
    #} end if
    mapping = dict()
    self.GetLine()
    #ExtremeDebugMsg(self, "Finished? %s" % getattr(self, "finished", False))
    #ExtremeDebugMsg(self, "CURR LINE: %s" % self.curr_line)
    while (not getattr(self, "finished", False) and ID_UPDATE_HEADER != self.curr_line): #{
      #ExtremeDebugMsg(self, "  Processing update: %s" % self.curr_line)
      (old_id, new_id) = self.curr_line.split(", ")[:2]
      if ("." in old_id): #{
        old_id = old_id.split(".")[0]
      #} end if
      if ("." in new_id): #{
        new_id = new_id.split(".")[0]
      #} end if
      #ExtremeDebugMsg(self, "    Old ID: %s, New ID: %s" % (old_id, new_id))
      if (old_id != new_id): #{
        if (old_id not in mapping): #{
          mapping[old_id] = set()
        #} end if
        mapping[old_id].add(new_id)
      #} end if
      self.GetLine()
    #} end while
    return mapping
  #} end def

  def GetLine(self): #{
    MAX_LOOPS = 100
    count = 0
    #self.GetRawLine()
    next_line = ""
    while (not getattr(self, "finished", False) and "" == next_line): #{
        #"" == getattr(self, "curr_line", "")): #{}
      #self.GetRawLine()
      next_line = self.stream.next()
      # next() should return at least "\n" for a blank line
      if ("" == next_line): #{
        #ExtremeDebugMsg(self, "STREAM COMPLETE")
        self.finished = True
        break
      #} end if
      next_line = CleanLine(next_line)
      #ExtremeDebugMsg(self, "UPDATE LINE: %s" % next_line)
      count += 1
      if (MAX_LOOPS < count): #{
        raise MainClassError("Too many loops while getting update line!")
      #} end if
      #ExtremeDebugMsg(self, "POLLING STREAM")
      # ret will be None until the process terminates
      #ret = self.stream.poll()
      #if (None != ret): #{
      #  ExtremeDebugMsg(self, "STREAM COMPLETE")
      #  self.finished = True
      #} end if
    #} end while
    self.curr_line = next_line
  #} end def

  #def GetRawLine(self): #{
  #  next_line = self.stream.next()
  #  # next() should return at least "\n" for a blank line
  #  if ("" == next_line): #{
  #    ExtremeDebugMsg(self, "STREAM COMPLETE")
  #    self.finished = True
  #  #} end if
  #  self.curr_line = CleanLine(next_line)
  #} end def
#} end class

#### EXCEPTION CLASSES ####
class CoverageParseError(MyError): #{
  pass
#} end class

class MainClassError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Use dwgsim to generate simulated reads for wildtype "
    "and event sequences.")
  #description_string = ("Associates each coverage value with a sequence "
  #    "and converts the coverage to an absolute number of reads for that "
  #    "sequence, then simulates reads from each sequence using the dwgsim "
  #    "tool.")
  args = [ "WT_ANNOT", "WT_SEQ", "WT_COV", "EVENT_SEQ", "EVENT_COV" ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  #parser.add_option("-A", "--cov-adjust",
  #    type="float", metavar="F",
  #    help="Add F to each coverage value read in from the coverages file "
  #      "(to ensure minimum coverage value) [default: %default]")
  parser.add_option("-R", "--read-length",
      type="int", metavar="N",
      help="Simulate reads Nnt long. [default: %default]")
  parser.add_option("-F", "--frag-length",
      type="int", metavar="N",
      help="Simulate read-pairs from fragments with a mean length of Nnt. "
        "[default: %default]")
  parser.add_option("--min-overlap",
      type="int", metavar="N",
      help="Ensure that for each event sequence, there is at least one read "
        "generated that overlaps the event split position by Nnt. "
        "[default: %default]")
  parser.add_option("--chr-filter",
      help="Only simulate reads from genes on chromosomes in the provided "
        "comma-delimited list.")
  dwgsim_group = OptionGroup(parser, "dwgsim Options")
  dwgsim_group.add_option("--err-rate",
      type="float", metavar="F",
      help="The per-base rate of sequencing errors. [default: %default]")
      #help="The per-base rate of sequencing errors, defaults to dwgsim "
      #  "default if unused.")
  dwgsim_group.add_option("--mut-rate",
      type="float", metavar="F",
      help="The rate of mutations. [default: %default]")
  dwgsim_group.add_option("--std-dev",
      type="float", metavar="F",
      help="The standard deviation of the fragment length, defaults to dwgsim "
        "default if unused.")
  dwgsim_group.add_option("--extra-dwgsim",
      help="A string containing any extra options to be used in the call to "
        "dwgsim.")
  parser.add_option_group(dwgsim_group)
  misc_group = OptionGroup(parser, "Miscellaneous Options")
  misc_group.add_option("--seed",
      type="float",
      help="The seed to use to initialize the random number generator. If no "
        "value is specified, current system time is used.")
  misc_group.add_option("--disable-profiling-timer",
                    action="store_true", dest="dpt",
                    help="Sometimes this script can hang when trying to spawn "
                         "child processes, due to the kernel's profiling "
                         "timer. Use this option to disable the profiling "
                         "timer if the script seems to be hanging.")
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
  parser.set_defaults(#cov_adjust=0.6,
                      read_length=75,
                      frag_length=200,
                      min_overlap=5,
                      err_rate=0.0037,
                      mut_rate=0,
                      dpt=False,
                      force=False,
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
  CheckFilePath(options.annot_path, "wildtype annotations file", path_errors)
  CheckFilePath(options.wseq_path, "wildtype sequences file", path_errors)
  CheckFilePath(options.wcov_path, "wildtype coverages file", path_errors)
  CheckFilePath(options.eseq_path, "event sequences file", path_errors)
  CheckFilePath(options.ecov_path, "event coverages file", path_errors)
  # get and check the output path
  #options.output_dir = GetOutDir(os.path.dirname(options.barnacle_path),
  #  "TASK_DESC")
  # create wildtype output paths
  (winput_dir, win_name) = os.path.split(options.wseq_path)
  options.woutput_dir = os.path.join(winput_dir, "simulated_reads")
  EnsureDirectoryExists(options.woutput_dir)
  wnreads_name = "%s.nreads_and_cov.txt" % os.path.splitext(win_name)[0]
  options.wnreads_path = os.path.join(options.woutput_dir, wnreads_name)
  wread1_seqs_name = "%s.r1.fq" % os.path.splitext(win_name)[0]
  options.wread1_seqs_path = os.path.join(options.woutput_dir, wread1_seqs_name)
  wread2_seqs_name = "%s.r2.fq" % os.path.splitext(win_name)[0]
  options.wread2_seqs_path = os.path.join(options.woutput_dir, wread2_seqs_name)
  # create event output paths
  (einput_dir, ein_name) = os.path.split(options.eseq_path)
  options.eoutput_dir = os.path.join(einput_dir, "simulated_reads")
  EnsureDirectoryExists(options.eoutput_dir)
  enreads_name = "%s.nreads_and_cov.txt" % os.path.splitext(ein_name)[0]
  options.enreads_path = os.path.join(options.eoutput_dir, enreads_name)
  eread1_seqs_name = "%s.r1.fq" % os.path.splitext(ein_name)[0]
  options.eread1_seqs_path = os.path.join(options.eoutput_dir, eread1_seqs_name)
  eread2_seqs_name = "%s.r2.fq" % os.path.splitext(ein_name)[0]
  options.eread2_seqs_path = os.path.join(options.eoutput_dir, eread2_seqs_name)
  if (not options.force): #{
    # check wildtype output paths
    CheckNewFilePath(options.wnreads_path, "wildtype read numbers file",
        path_errors)
    CheckNewFilePath(options.wread1_seqs_path,
        "wildtype read 1 sequences file", path_errors)
    CheckNewFilePath(options.wread2_seqs_path,
        "wildtype read 2 sequences file", path_errors)
    # check wildtype output paths
    CheckNewFilePath(options.enreads_path, "event read numbers file",
        path_errors)
    CheckNewFilePath(options.eread1_seqs_path, "event read 1 sequences file",
        path_errors)
    CheckNewFilePath(options.eread2_seqs_path, "event read 2 sequences file",
        path_errors)
  #} end if
  if (opts_good and 0 == len(path_errors)): #{
    #CheckDirPath(options.output_dir, "output", path_errors,
    #  create=True, replace=options.force)
    # get the log file name
    options.log_file_name = GetLogPath(options.wseq_path,
      "simulate_reads", options.woutput_dir)
  #} end if
  if (0 < len(path_errors)): #{
    ErrMsg("Errors in input arguments:\n  %s" % "\n  ".join(path_errors))
  #} end if
  # seed the random number generator
  if (None == options.seed): #{
    options.seed = IntFloor(time.time())
  #} end if
  random.seed(options.seed)
  # set up the chromosome filter list, if one has been given
  if (None != options.chr_filter): #{
    options.chr_filter = set(map(NormalizeChrID,
      options.chr_filter.split(",")))
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
    options.annot_path = EnsureAbsPath(args[0])
    options.wseq_path = EnsureAbsPath(args[1])
    options.wcov_path = EnsureAbsPath(args[2])
    options.eseq_path = EnsureAbsPath(args[3])
    options.ecov_path = EnsureAbsPath(args[4])
    if (CheckPaths(options)): #{
      try: #{
        main_class_object = MainClassCls(options)
        WriteCommand(main_class_object, sys.argv)
        main_class_object.Run()
      except (MyError), e:
        ErrMsg("ERROR while simulating reads from wildtype and "
            "event sequences:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify the path to a fasta sequences file "
      "(SEQ_FILE); and the path to a file listing coverages (COV_FILE).")
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
