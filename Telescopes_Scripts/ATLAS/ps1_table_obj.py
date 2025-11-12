
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import FITSFixedWarning

import random
import string
import time

import mastcasjobs
from AstroColour.AstroColour import RGB

import warnings # To ignore our problems
warnings.filterwarnings('ignore', category=FITSFixedWarning)

user = "syndiff"
pwd = "PS1SynDiff"

def ps1_casjobs(ra, dec, radius, maglim=20, quick = False):
    
    name = ''.join(random.choices(string.ascii_letters, k=8))
    
    query = f"""
            SELECT o.objID,
            o.raMean, o.decMean,
            o.qualityFlag,
            o.gMeanPSFMag, o.gMeanPSFMagErr, o.gMeanPSFMagNpt,
            o.rMeanPSFMag, o.rMeanPSFMagErr, o.rMeanPSFMagNpt,
            o.iMeanPSFMag, o.iMeanPSFMagErr, o.iMeanPSFMagNpt,
            o.zMeanPSFMag, o.zMeanPSFMagErr, o.zMeanPSFMagNpt,
            o.yMeanPSFMag, o.yMeanPSFMagErr, o.yMeanPSFMagNpt,
            o.rMeanKronMag, o.rMeanKronMagErr,
            o.iMeanKronMag, o.iMeanKronMagErr,
            o.nDetections,
            o.gFlags, o.gQfPerfect,
            o.rFlags, o.rQfPerfect,
            o.iFlags, o.iQfPerfect,
            o.zFlags, o.zQfPerfect,
            o.yFlags, o.yQfPerfect
            INTO mydb.[{name}]
            from fGetNearbyObjEq({ra}, {dec}, {radius}*60) as x
            JOIN MeanObjectView o on o.ObjID=x.ObjId
            LEFT JOIN StackObjectAttributes AS soa ON soa.objID = x.objID
            WHERE o.nDetections>5
            AND soa.primaryDetection>0
            AND o.iMeanPSFMag < 20
            """

    jobs = mastcasjobs.MastCasJobs(username=user, password=pwd, context="PanSTARRS_DR2")
    
    if quick:
        table = quick_query(jobs, query)
        if table is None:
            print('Error in casjobs quick query, doing solid query')
            table = solid_query(jobs, query, name)
    else:
        table = solid_query(jobs, query, name)
        
    return table

def quick_query(jobs, query):
    try:
        results = jobs.quick(query, task_name="python cone search")
        return results
        # return results.to_pandas()
    except:
        return None

def solid_query(jobs, query, name):

    results = jobs.submit(query, task_name="python cone search")
    job_status = jobs.status(results)
    
    stsus = False
    print("Finding...", end="", flush=True)
    while stsus == False:
        time.sleep(5)
        print(".", end="", flush=True)
        job_status = jobs.status(results)
        if job_status[1].lower() == 'finished':
            stsus = True
        elif job_status[1].lower() == 'failed':
            print('Error in casjobs query')
            stsus = None
        else:
            stsus = False
            
    if stsus is True:
        
        print("Query finished successfully. Retrieving table...")

        table = jobs.get_table(f"mydb.{name}")
        return table.to_pandas()
    else:
        print('Error in casjobs query')
        
