#! /usr/bin/env python
"""
create_gene_feature_coords.py

Created by Lucas Swanson
Copyright (c) 2012 Canada's Michael Smith Genome Sciences Centre. All rights reserved.
"""

# import standard modules
from optparse import OptionParser, OptionGroup
import os, re, sys, time, traceback

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
  NormalizeChrID, NonStandardChr)
from utils.messages import ErrMsg, LogMsg, DebugMsg, ExtremeDebugMsg
from utils.files_paths import (CheckFilePath, CheckDirPath, EnsureAbsPath,
  GetOutDir, FileBoxCls)
from parsers.genes.annotation import (GeneAnnotationParserCls, ANNOT_TYPES)

# CONSTANTS
ES_SUCCESS = 0
ES_RUN_ERR = 1
ES_OPT_ERR = 2
ES_PATH_ERR = 3
ES_EXCEPTION = 4
MSG_SUCCESS = "GENE FEATURE SUCCESS"
MSG_FAIL = "GENE FEATURE FAIL"
EXON = "exon"
INTRON = "intron"
UTR = "utr"

class ExonCoordsCreatorCls: #{
  def __init__(self, options): #{
    SetupMainClass(self, options)
  #} end def

  def __del__(self): #{
    CloseLogFile(self)
  #} end def

  def CreateExonCoordsFile(self): #{
    LogMsg(self, "Creating gene feature coordinates file...\n%s" %
      self.options.output_path)
    start = time.time()
    # get the annotations type
    #if (None == self.options.annotations_type): #{
    #  annots_type = GetAnnotationsType(self.options.input_file)
    #else:
    #  annots_type = self.options.annotations_type
    #} end if
    #parse_annot = GetAnnotationsParseFunction(annots_type, self.log_info)
    # open the input and output files
    #in_file = FileBoxCls(self.options.input_file, "r",
    #  "cannot open gene annotations file")
    # get the annotations parser
    parser = GeneAnnotationParserCls(self.options.input_file,
      type=self.options.annotations_type, log_info=self.log_info)
    out_file = FileBoxCls(self.options.output_path, "w",
      "cannot create exon coordinates file")
    # iterate through the gene annotations file
    LogMsg(self, "Parsing gene annotations file...")
    #for line in in_file: #{}
    num_genes = 0
    for annot in parser: #{
      #DebugMsg(self, line)
      #annot = FixAnnotation(parse_annot(line))
      annot = FixAnnotation(annot)
      if (self.SkipAnnot(annot)): #{
        continue
      #} end if
      # check whether the gene is non-coding
      if (annot.non_coding): #{
        DebugMsg(self, "Non-coding transcript %s %i-%i" %
          (annot.name, annot.cdsStart, annot.cdsEnd))
      #} end if
      self.WriteFeatures(annot, out_file)
      num_genes += 1
    #} end for
    out_file.Close()
    DebugMsg(self, "Wrote feature coordinates for %i genes." % num_genes)
    LogMsg(self, "Time spent creating exon coordinates file: %s" %
      TimeSpent(start))
  #} end def

  def SkipAnnot(self, annot): #{
    # filter chromosomes, if the option is being used
    if (self.options.filter_chromosomes): #{
      #chr_patt = r"\A(chr)?([1-9][0-9]*|[XY]|MT?)\Z"
      #if (None == re.search(chr_patt, annot.chrom)): #{}
      if (NonStandardChr(annot.chrom)): #{
        DebugMsg(self, "Skipping exon in chromosome %s" % annot.chrom)
        return True
      #} end if
    #} end if
    return False
  #} end def

  def WriteFeatures(self, annot, out_file): #{
    prev = None
    # iterate through the exons for the gene
    #num_exons = len(annot.exons)
    exon_list = annot.exons
    if ("-" == annot.strand and 1 < len(annot.exons) and
        annot.exons[0][0] < annot.exons[1][0]): #{
      exon_list = reversed(annot.exons)
    #} end if
    for i, coords in enumerate(exon_list): #{
      if (None == coords): #{
        DebugMsg(self, "Skipping empty exon %i" % i)
        continue
      #} end if
      #(e_start, e_end) = exon
      exon = OutputFeatureCls(annot, EXON, coords, i)
      msg = "Exon %i start: %i, end: %i " % (i, exon.start, exon.end)
      if (None == prev): #{
        msg += "prev: None"
      else:
        msg += "prev_s: %i, prev_e: %i" % (prev[0], prev[1])
      #} end if
      DebugMsg(self, msg)
      feature_strings = list()
      if (self.options.include_introns and None != prev): #{
        #intron_data = list([
        #  annot.chrom,
        #  "%i" % (prev_end+1),
        #  "%i" % (exon.start-1),
        #  "%s.intron-%i" % (annot.name, i-1),
        #])
        #feature_strings.append(" ".join(intron_data))
        #DebugMsg(self, "  intron: %s" % " ".join(intron_data))
        if ("+" == annot.strand): #{
          intron = OutputFeatureCls(annot, INTRON,
            (prev[1]+1,exon.start-1), i-1)
        else:
          intron = OutputFeatureCls(annot, INTRON,
            (exon.end+1,prev[0]-1), i-1)
        #} end if
        feature_strings.append(intron.ToString())
        DebugMsg(self, "  intron: %s" % intron.ToString())
      #} end if
      prev = coords
      # if the gene is non-coding, just write the exon
      if (annot.non_coding): #{
        #exon_data = list([
        #  annot.chrom,
        #  "%i" % exon.start,
        #  "%i" % exon.end,
        #  "%s-%i" % (annot.name, i),
        #])
        #feature_strings.append(" ".join(exon_data))
        #DebugMsg(self, "  non-coding: %s" % " ".join(exon_data))
        feature_strings.append(exon.ToString())
        DebugMsg(self, "  non-coding: %s" % exon.ToString())
      # if the exon ends before the CDS start or
      # the exon starts after the CDS end,
      # the full exon is a UTR
      elif (exon.end < annot.cdsStart or annot.cdsEnd < exon.start):
        #exon_data = list([
        #  annot.chrom,
        #  "%i" % exon.start,
        #  "%i" % exon.end,
        #  "%s.UTR-%i" % (annot.name, i),
        #])
        #feature_strings.append(" ".join(exon_data))
        #DebugMsg(self, "  full UTR: %s" % " ".join(exon_data))
        exon.type = UTR
        feature_strings.append(exon.ToString())
        DebugMsg(self, "  full UTR: %s" % exon.ToString())
      else:
        # if the exon starts before the CDS start and
        # ends after the CDS start,
        # the first part of the exon is a UTR
        if (exon.start < annot.cdsStart): #{
          #utr_data = list([
          #  annot.chrom,
          #  "%i" % exon.start,
          #  "%i" % (annot.cdsStart-1),
          #  "%s.UTR-%i" % (annot.name, i),
          #])
          #feature_strings.append(" ".join(utr_data))
          utr = OutputFeatureCls(annot, UTR,
            (exon.start, annot.cdsStart-1), i)
          feature_strings.append(utr.ToString())
          exon.start = annot.cdsStart
          #DebugMsg(self, "  UTR start: %s\n  New start: %i" %
          #  (" ".join(utr_data), exon.start))
          DebugMsg(self, "  UTR start: %s\n  New start: %i" %
            (utr.ToString(), exon.start))
        #} end if
        # if the exon starts before the CDS end and
        # ends after the CDS end,
        # the second part of the exon is a UTR
        if (annot.cdsEnd < exon.end): #{
          #exon_data = list([
          #  annot.chrom,
          #  "%i" % exon.start,
          #  "%i" % annot.cdsEnd,
          #  "%s-%i" % (annot.name, i),
          #])
          #feature_strings.append(" ".join(exon_data))
          utr = OutputFeatureCls(annot, UTR, (annot.cdsEnd+1, exon.end), i)
          exon.end = annot.cdsEnd
          feature_strings.append(exon.ToString())
          #utr_data = list([
          #  annot.chrom,
          #  "%i" % (annot.cdsEnd+1),
          #  "%i" % exon.end,
          #  "%s.UTR-%i" % (annot.name, i),
          #])
          #feature_strings.append(" ".join(utr_data))
          #DebugMsg(self, "  exon start: %s\n  UTR end: %s" %
          #  (" ".join(exon_data), " ".join(utr_data)))
          feature_strings.append(utr.ToString())
          DebugMsg(self, "  exon start: %s\n  UTR end: %s" %
            (exon.ToString(), utr.ToString()))
        # if the exon starts after the CDS start and
        # ends before the CDS end,
        # the full exon is really an exon
        elif (exon.start <= exon.end):
          #exon_data = list([
          #  annot.chrom,
          #  "%i" % exon.start,
          #  "%i" % exon.end,
          #  "%s-%i" % (annot.name, i),
          #])
          #feature_strings.append(" ".join(exon_data))
          #DebugMsg(self, "  full exon: %s" % " ".join(exon_data))
          feature_strings.append(exon.ToString())
          DebugMsg(self, "  full exon: %s" % exon.ToString())
        else:
          raise ExonCoordsError("cannot determine exon type: "
            "%s: CDS:%i-%i, Exon-%i:%i-%i" % (annot.name,
            annot.cdsStart, annot.cdsEnd, i, exon.start, exon.end))
        #} end if
      #} end if
      out_file.WriteLine("%s" % "\n".join(feature_strings))
    #} end for
  #} end def
#} end class

class OutputFeatureCls: #{
  def __init__(self, annot, type, coords, index): #{
    self.transcript_name = annot.name
    self.non_coding = annot.non_coding
    if (hasattr(annot, "alias")): #{
      self.gene_name = annot.alias
    else:
      self.gene_name = None
    #} end if
    self.type = type
    self.index = index
    self.chrom = annot.chrom
    (self.start, self.end) = coords
  #} end def

  def ToString(self): #{
    id = "%s" % self.transcript_name
    if (self.non_coding): #{
      id += ".non_coding"
    #} end if
    if (None != self.gene_name): #{
      id += ">%s" % self.gene_name
    #} end if
    id += "~%s~%i" % (self.type, self.index)
    return " ".join([self.chrom, "%i" % self.start, "%i" % self.end, id])
  #} end def
#} end class

def FixAnnotation(annot, use_chr=True): #{
  if (None == annot): #{
    return None
  #} end if
  # standardize the chromosome name
  #annot.chrom = annot.chrom.replace("MT","M")
  #if (not annot.chrom.startswith("chr")): #{
  #  annot.chrom = "chr%s" % annot.chrom
  #} end if
  annot.chrom = NormalizeChrID(annot.chrom, use_chr=use_chr)
  # replace spaces and tildes in the name and alias
  annot.name  = annot.name.replace(" ", "_").replace("~","_")
  # some annotation files do not have an "alias" field (e.g. "ensg")
  if (hasattr(annot, "alias") and None != annot.alias): #{
    annot.alias = annot.alias.replace(" ", "_").replace("~","_")
  #} end if
  # note that start positions need to be adjusted
  # (but exon start positions were adjusted during parsing)
  annot.txStart = int(annot.txStart)+1
  annot.txEnd = int(annot.txEnd)
  if (None == annot.cdsStart): #{
    annot.cdsStart = annot.txStart
  else:
    annot.cdsStart = int(annot.cdsStart)+1
  #} end if
  if (None == annot.cdsEnd): #{
    annot.cdsEnd = annot.txEnd
  else:
    annot.cdsEnd = int(annot.cdsEnd)
  #} end if
  if (hasattr(annot, "exonCount") and None != annot.exonCount): #{
    annot.exonCount = int(annot.exonCount)
  #} end if
  # check whether the gene is non-coding
  if (annot.cdsStart >= annot.cdsEnd): #{
    annot.non_coding = True
  else:
    annot.non_coding = False
  #} end if
  return annot
#} end def

#### EXCEPTION CLASSES ####
class ExonCoordsError(MyError): #{
  pass
#} end class

def SetupOptionsParser(): #{
  description_string = ("Extract exon coordinates from a gene annotation "
    "file and create a file holding the exon coordinates formatted for the "
    "overlap code, with non-coding genes and UTRs marked.")
  args = [ "INPUT_FILE", "OUTPUT_DIR", ]
  usage_string = "%prog " + " ".join(args) + " [ OPTIONS ]"
  parser = OptionParser(description=description_string,
                        version="%prog " + VERSION,
                        usage=usage_string)
  parser.num_args = len(args)
  parser.add_option("--annotations-type",
                    help="The type of annotations file being used for the "
                         "gene names, e.g.: %s. " % ", ".join(ANNOT_TYPES) +
                         "The code will try to automatically determine the "
                         "annotations type if this option is not used.")
  parser.add_option("--filter-chromosomes",
                    action="store_true",
                    help="Only write out exons for genes on standard "
                         "chromosomes: chr<I>, chrX, chrY, chrM (where "
                         "<I> is any integer). [default]")
  parser.add_option("--all-chromosomes",
                    action="store_false", dest="filter_chromosomes",
                    help="Write out exons for genes on all chromosomes.")
  parser.add_option("--no-introns",
                    action="store_false", dest="include_introns",
                    help="Write only exon coordinates, not introns.")
  parser.add_option("--introns",
                    action="store_true", dest="include_introns",
                    help="Write both exon and intron coordinates. [default]")
  parser.add_option("-f", "--force",
                    action="store_true",
                    help="Force filtering to take place, overwriting the exon "
                         "coordinates file if it already exists.")
  parser.add_option("-d", "--debug",
                    action="store_true",
                    help="Print debug information while the program runs.")
  parser.add_option("--extreme-debug",
                    action="store_true", dest="extreme_debug",
                    help="Print extremely in-depth debug information while "
                      "the program runs. Not recommended for large jobs.")
  parser.set_defaults(filter_chromosomes=True,
                      include_introns=True,
                      force=False,
                      debug=False,
                      extreme_debug=False)
  return parser
#} end def

def CheckPaths(options): #{
  opts_good = True
  if (None != options.annotations_type and
      options.annotations_type.lower() not in ANNOT_TYPES):
    ErrMsg("\"%s\" is not a valid annotations type. Valid types are: %s" %
            ", ".join(ANNOT_TYPES))
    opts_good = False
  #} end if
  path_errors = list()
  CheckFilePath(options.input_file, "gene annotations", path_errors)
  if (opts_good and 0 == len(path_errors)): #{
    CheckDirPath(options.output_dir, "output", path_errors, create=True)
  #} end if
  if (opts_good and 0 == len(path_errors)): #{
    input_file_name = os.path.basename(options.input_file)
    (input_root, ext) = os.path.splitext(input_file_name)
    output_root = "%s.exons" % input_root
    if (options.include_introns): #{
      output_root += ".introns"
    #} end if
    if (options.filter_chromosomes): #{
      output_root += ".std_chr"
    #} end if
    output_file_name = output_root + ".bed"
    options.output_path = os.path.join(options.output_dir, output_file_name)
    if (options.debug): #{
      ErrMsg("Output path: %s" % options.output_path)
    #} end if
    if (not options.force and os.path.exists(options.output_path)): #{
      path_errors.append("Exon coordinates file already exists, please "
        "remove it or use --force option to overwrite it: \"%s\"." %
        options.output_path)
    #} end if
    # set the log-file name
    options.log_file_name = GetLogPath(options.input_file,
      "exon_coords", options.output_dir)
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
    options.input_file = EnsureAbsPath(args[0])
    options.output_dir = EnsureAbsPath(args[1])
    if (CheckPaths(options)): #{
      try: #{
        exon_coords_creator = ExonCoordsCreatorCls(options)
        WriteCommand(exon_coords_creator, sys.argv)
        exon_coords_creator.CreateExonCoordsFile()
      except (MyError), e:
        ErrMsg("ERROR while creating exon coordinates file:\n  %s" % e)
        return ES_RUN_ERR
      #} end try
    else:
      return ES_PATH_ERR
    #} end if
  else:
    parser.error("you must specify the path to a GTF or UCSC gene annotations "
      "file (INPUT_FILE); and the path to a directory to hold the output "
      "(OUTPUT_DIR).")
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
