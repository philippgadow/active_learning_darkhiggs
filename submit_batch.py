import os
import re
import sys
import argparse
import logging
from pprint import pformat, pprint
from batch_management.condor_handler import CondorHandler


"""
Script for submitting jobs to a HTCondor batch system.
"""


def getArgumentParser():
  parser = argparse.ArgumentParser()
  parser.add_argument("--batch_dir", default=None, help="Directory for batch submission files")
  parser.add_argument("--out_dir", default=None, help="Output directory for results")
  parser.add_argument('--grid', default='grid.csv')
  parser.add_argument("--msglevel", default="info", choices=["info", "debug", "error"], help="Message output detail level.")
  parser.add_argument("--start_dsid", default=100000, help="Start to incremental DSID generation")
  return parser


def ensureDirExists(item):
  """Create a directory if it does not already exist."""
  if not os.path.exists(item):
    os.makedirs(item)


def ensureDirsExist(dir_list):
  """Create a list of directories if it does not already exist."""
  for item in dir_list:
    ensureDirExists(item)


def setupCondorHandler(batch_dir) :
  """Setup HTCondor handle to submit jobs."""
  logging.debug("Setting up HTCondor handler (only done once)...")
  batch_path = os.path.join(batch_dir, "batch")
  log_path = os.path.join(batch_dir, "batch_logs")
  ensureDirsExist([batch_path, log_path])

  handler = CondorHandler(batch_path, log_path)
  handler['runtime'] = '10800' # in seconds, 10800 = 3h
  handler['memory'] = '2GB'
  handler['cpu'] = '1'
  handler['project'] = 'af-atlas'  # for DESY NAF
  # handler['jobflavour'] = 'workday'  # only on lxplus
  # handler['container'] = ''
  return handler


def submitJob(batch_dir, dsid, mzp, mdh, mdm, gq, gx):
  """Submit job using htcondor python bindings."""
  logging.debug("Submitting job: {0}".format(dsid))

  workdir = os.getcwd()
  tag = "job_{0}".format(dsid)

  if not hasattr(submitJob, "handler"):
      submitJob.handler = setupCondorHandler(batch_dir)

  # command on htcondor host:
  command = "cd {workdir} && bash {workdir}/run_workflow.sh {dsid} {mzp} {mdh} {mdm} {gq} {gx}".format(
    workdir=workdir, dsid=dsid, mzp=mzp, mdh=mdh, mdm=mdm, gq=gq, gx=gx
    )

  # submitJob.handler.activate_testmode()
  submitJob.handler.send_job(command, tag)


def main():
  # get arguments from command line
  args = getArgumentParser().parse_args()

  # configure logging output level
  output_level = {
    "info": logging.INFO,
    "debug": logging.DEBUG ,
    "error": logging.DEBUG
  }
  logging.basicConfig(level=output_level[args.msglevel])

  # parameters (eventually to be outsourced to a config file)
  # dsid, mzp,mdh,mdm,gq,gx
  job_parameters = []
  with open(args.grid) as f:
    for l in f: job_parameters.append(l.strip())
  pprint(job_parameters)

  # create directories for condor submission and output
  logging.info("Preparing batch job submission...")
  batch_dir = args.batch_dir if args.batch_dir else os.path.join(os.getcwd(), "condor")
  ensureDirsExist([batch_dir])
  
  dsid = args.start_dsid
  # submit jobs
  logging.info("Submitting jobs...")
  for data in job_parameters:
    if not data.split(): continue
    _, mzp, mdh, mdm, gq, gx = data.split(',')
    mzp = round(float(mzp))
    mdh = round(float(mdh))
    mdm = round(float(mdm))
    if not os.path.isfile(os.path.join('data', str(dsid), 'histograms.root')):
      submitJob(batch_dir, str(dsid), mzp, mdh, mdm, gq, gx)
    dsid += 1

if __name__ == '__main__':
  main()
