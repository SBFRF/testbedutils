import os, logging
import datetime as DT


def makeCMTBfileStructure(path_prefix, date_str):
    """checks and makes sure there is a folder structure that can beused for file storage"""
    if not os.path.exists(os.path.join(path_prefix, date_str)):  # if it doesn't exist
        os.makedirs(os.path.join(path_prefix, date_str))  # make the directory
    if not os.path.exists(os.path.join(path_prefix, date_str, 'figures')):
        os.makedirs(os.path.join(path_prefix, date_str, "figures"))
    print('use function for making files ')
    print("OPERATIONAL files will be place in {0} folder".format(os.path.join(path_prefix, date_str)))

def logFileLogic(outDataBase, version_prefix, startTime, endTime, log=True):
    """Checks and makes log file

    Args:
        outDataBase: string for working directory
        version_prefix: version prefix for model
        startTime: start time of project simulation
        endTime: end time of project simulation
        log(bool): turn logs on or not (default=True)

    Returns:

    """


    LOG_FILENAME = outDataBase + 'logs/cmtb_BatchRun_Log_%s_%s_%s.log' % (version_prefix, startTime, endTime)
    if log is True:
        try:
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
        except IOError:
            os.makedirs(outDataBase+'logs')
            logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)
        logging.debug('\n-------------------\nTraceback Error Log for:\n\nSimulation Started: %s\n-------------------\n'
                  % (DT.datetime.now()))
    # ____________________________________________________________
    return LOG_FILENAME
def displayStartInfo(projectStart, projectEnd, version_prefix, LOG_FILENAME, model):
    print('\n-\n-\nMASTER WorkFLOW for {} SIMULATIONS\n-\n-\n'.format(model))
    print('Batch Process Start: %s     Finish: %s '% (projectStart, projectEnd))
    print('The batch simulation is Run in %s Version' % version_prefix)
    print('Check for simulation errors here %s' % LOG_FILENAME)

    print('------------------------------------\n\n************************************\n\n------------------------------------\n\n')

def checkVersionPrefix(model, version_prefix):
    """ a function to check if model prefix is already programmed into structure

    Args:
        model: a model name string
        version_prefix: version prefix string

    Returns:
        None

    """
    cmsStrings = ['base']
    ww3Strings = ['base']
    stwaveStrings= ['HP', 'FP', 'CB', 'CBKF']

    if model.lower() in ['cms']:
        modelList = cmsStrings
    elif model.lower() in ['ww3']:
        modelList = ww3Strings
    elif model.lower() in ['stwave']:
        modelList = stwaveStrings
    else:
        raise NotImplementedError, 'Check model is programmed'
    checkString = 'Your model is not in version Prefix list {}'.format(modelList)
    assert version_prefix.lower() in cmsStrings, checkString