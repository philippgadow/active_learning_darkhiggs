import os
import re
import sys
import argparse
import logging
from pprint import pformat
from batch_management.condor_handler import CondorHandler


"""
Script for submitting jobs to a HTCondor batch system.
"""


def getArgumentParser():
  parser = argparse.ArgumentParser()
  parser.add_argument("--batch_dir", default=None, help="Directory for batch submission files")
  parser.add_argument("--out_dir", default=None, help="Output directory for results")
  parser.add_argument("--msglevel", default="info", choices=["info", "debug", "error"], help="Message output detail level.")
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
  handler['runtime'] = '7200' # in seconds, 7200s = 2h
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
  # write output to /tmp dir on batch machine and not directly to output folder
  # to enhance performance by avoiding long access of shared filesystems
  # command = "rm -rf /tmp/{0} && mkdir -p /tmp/{0} && ".format(tag)
  # command = "export DSID={0} && ".format(dsid)
  # command = "export MZP={0} && ".format(mzp)
  # command = "export MDH={0} && ".format(mdh)
  # command = "export MDM={0} && ".format(mdm)
  # command = "export GQ={0} && ".format(gq)
  # command = "export GX={0} && ".format(gx)
  command = "cd {workdir} && bash {workdir}/run_workflow.sh {dsid} {mzp} {mdh} {mdm} {gq} {gx}".format(
    workdir=workdir, dsid=dsid, mzp=mzp, mdh=mdh, mdm=mdm, gq=gq, gx=gx
    )
  # command += "cp /tmp/{0}/* {1} && rm -rf /tmp/{0}".format(tag, out_dir)

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
  job_parameters = [
    "100000,500,50,200,0.25,1.0",
    "100001,500,70,200,0.25,1.0",
    "100002,500,90,200,0.25,1.0",
    "100003,500,110,200,0.25,1.0",
    "100004,500,130,200,0.25,1.0",
    "100005,500,150,200,0.25,1.0",
    "100006,1000,50,200,0.25,1.0",
    "100007,1000,70,200,0.25,1.0",
    "100008,1000,90,200,0.25,1.0",
    "100009,1000,110,200,0.25,1.0",
    "100009,1000,130,200,0.25,1.0",
    "100010,1000,150,200,0.25,1.0",
    "100011,1500,50,200,0.25,1.0",
    "100012,1500,70,200,0.25,1.0",
    "100013,1500,90,200,0.25,1.0",
    "100014,1500,110,200,0.25,1.0",
    "100015,1500,130,200,0.25,1.0",
    "100016,1500,150,200,0.25,1.0",
    "100017,2000,50,200,0.25,1.0",
    "100018,2000,70,200,0.25,1.0",
    "100019,2000,90,200,0.25,1.0",
    "100020,2000,110,200,0.25,1.0",
    "100021,2000,130,200,0.25,1.0",
    "100022,2000,150,200,0.25,1.0",
    "100023,2500,50,200,0.25,1.0",
    "100024,2500,70,200,0.25,1.0",
    "100025,2500,90,200,0.25,1.0",
    "100026,2500,110,200,0.25,1.0",
    "100027,2500,130,200,0.25,1.0",
    "100028,2500,150,200,0.25,1.0",
    "100029,3000,50,200,0.25,1.0",
    "100030,3000,70,200,0.25,1.0",
    "100031,3000,90,200,0.25,1.0",
    "100032,3000,110,200,0.25,1.0",
    "100033,3000,130,200,0.25,1.0",
    "100034,3000,150,200,0.25,1.0",
    "100035,3500,50,200,0.25,1.0",
    "100036,3500,70,200,0.25,1.0",
    "100037,3500,90,200,0.25,1.0",
    "100038,3500,110,200,0.25,1.0",
    "100039,3500,130,200,0.25,1.0",
    "100040,3500,150,200,0.25,1.0"
  ]

  # create directories for condor submission and output
  logging.info("Preparing batch job submission...")
  batch_dir = args.batch_dir if args.batch_dir else os.path.join(os.getcwd(), "condor")
  ensureDirsExist([batch_dir])

  # submit jobs
  logging.info("Submitting jobs...")
  for data in job_parameters:
    dsid, mzp, mdh, mdm, gq, gx = data.split(',')
    submitJob(batch_dir, dsid, mzp, mdh, mdm, gq, gx)


if __name__ == '__main__':
  main()
