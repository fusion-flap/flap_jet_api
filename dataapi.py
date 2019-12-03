# -*- coding: utf-8 -*-
"""

.. module:: FLAP/***
   :platform: UNIX
   :synopsis: Signal reading for JET

.. moduleauthor:: M. Vecsei


"""
import ppf
import getdat as jpf
import os
import pwd
import flap
import numpy as np
from datetime import date

def getsignal(exp_id, source, no_data = False, options={}):
    """ Signal reading function for the JET API
    Does not require using the FLAP storage
    Input: same as get_data()
    Output: DataObject for the given source
    """
    # Function used for getting either ppf or jpf data from the server
    # 
    options_default = {"Sequence": 0,
                       "UID": "jetppf",
                       "fix0": 0,
                       "reshape": 0,
                       "no_x": 0,
                       "no_t": 0,
                       "Only Info": 0,
                       "nwds": 0,
                       "Check Time Equidistant": False,
                       "Cache Data": False}
    options = {**options_default, **options}
    
    #checking if the data is already in the cache
    curr_path = os.path.dirname(os.path.abspath(__file__))
    location = os.path.sep.join(curr_path.split(os.path.sep))
    split_source = source.split("/")
    filename = os.path.join(os.path.sep.join([location,"cached"]),
                            ("_".join(split_source)+"-"+options["UID"]+\
                            "-"+str(exp_id)+".hdf5").lower())
    print(filename)
    if os.path.exists(filename):
        print("ok")
        return flap.load(filename)

    # obtaining the data from the server
    if split_source[0].upper() == "PPF":
        raw_ppf = \
                ppf.ppfdata(exp_id,split_source[2],split_source[1],
                            seq=options["Sequence"], uid=options["UID"],
                            fix0=options["fix0"], reshape=options["reshape"],
                            no_x=options["no_x"], no_t=options["no_t"],
                            no_data=options["Only Info"])
        [data,x,t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier] = raw_ppf
        if ier != 0:
            raise ValueError("Error reading "+source+" for shot "+str(exp_id)+
                             ", UID "+options["UID"]+":\nIER "+str(ier)+": "+desc)
    elif split_source[0].upper() == "JPF":
        raw_jpf = jpf.getdat("/".join(split_source[1:]),
                             exp_id, nwds = options["nwds"])
        [data, t, nwds, desc, dunits, ier] = raw_jpf
        options["no_x"] = 0
        x=[]
        tunits = ""
        if ier != 0:
            raise ValueError("Error reading "+source+" for shot "+str(exp_id)+
                             ":\nIER "+str(ier))
    else:
        raise ValueError("Source should be either PPF or JPF")
    
    # converting data to flap DataObject
    coordinates = []
    info = "Obtained at "+str(date.today())+"\n"
    if options["no_t"] == 0 and (len(t) != 0):
        tunit = ["secs", "s", "seconds"]
        if tunits.lower() in tunit:
            tunits = "Second"
        equidist = False
        if options['Check Time Equidistant'] is True:
            timesteps = np.array(t[1:])-np.array(t[0:-1])
            equidistance = np.linalg.norm(timesteps-timesteps[0])/np.linalg.norm(timesteps)
            if equidistance < 1e-6:
                info = info + "Time variable is taken as equidistant to an accuracy of "+\
                       str(equidistance)+"\n"
                time_coord = flap.Coordinate(name='Time', unit=tunits, start = t[0],
                                             shape=t.shape, step=np.mean(timesteps),
                                             mode=flap.CoordinateMode(equidistant=True),
                                             dimension_list=[0])
                equidist = True
        if equidist is False:
            if options['Check Time Equidistant'] is True:
                info = info+"Time variable is not equidistant, deviation: "+str(equidistance)+"\n"
            time_coord = flap.Coordinate(name='Time', unit=tunits, values=t,
                                         shape=t.shape,
                                         mode=flap.CoordinateMode(equidistant=False),
                                         dimension_list=[0])
        coordinates.append(time_coord)
    if options["no_x"] == 0 and (len(x) != 0):
        x_coord = flap.Coordinate(name='X', unit=xunits, values=x,
                                     shape=x.shape,
                                     mode=flap.CoordinateMode(equidistant=False),
                                     dimension_list=[0])
        coordinates.append(x_coord)
    if options["Only Info"] == 0 and (len(data) != 0):
        unit = flap.Unit(name=source,unit=dunits)
        signal = flap.DataObject(data_array=data, data_unit=unit,
                                     coordinates=coordinates, exp_id=str(exp_id),
                                     data_title=desc, data_shape=data.shape,
                                     info=info)
        if "comm" in locals():
            signal.info = signal.info+comm+"\n"
        if no_data is True:
            signal.data = None
            
        if options["Cache Data"] is True:
            flap.save(signal, filename)
        return signal
    else:
        # if one is only ineterested in the metadata
        return raw_ppf

def get_data(exp_id=None, data_name=None, no_data=False, options={}, coordinates=None, data_source=None):
    """ Signal reading function for the JET API
    data_name: should be a string, e.g "PPF/ZPRO/EFIT" or "JPF/DH/Y6-EMITER<VLT"
    exp_id: Experiment ID, '95531'
    no_data: If set to True, than only the coordinates and the metadata of the
             signal is returned
    options:
        Time Start and Time Stop times can be given (see exp_id) e.g.
         options={
         Sequence - the sequence number, if undefined, than the latest data for a pulse is obtained
         UID - the user id under which the ppf is to be found. By default set to jetppf
               if the value is set to "curr_user", than the UID is set to the current user's ID
         no_x, no_t, Only Info - If no_x and/or no_t is set, then the corresponding
                               coordinate is not returned. If "Only Info" is
                               set then both x, t and the data is not returned,
                               only the metadata.
        Check Time Equidistant - if set and True, then the data is checked if it is equidistant in time.
                              it is assumed to be equidistant, if the difference between the timesteps relatively small
        nwds, fix0, reshape - Check jet ppf guide for info about these
        Cache Data - if True, the data will be cached in a file to obtain later
    Output: DataObject for the data_name and exp_id
    """
    options_default = {
            "Sequence": 0,
            "UID": "jetppf",
            "no_x": 0,
            "no_t": 0,
            "Only Info": 0,
            "Check Time Equidistant": False,
            "Cache Data": False,
            "fix0": 0,
            "reshape": 0,
            "nwds": 0}
    options = {**options_default, **options}
    
    if options["UID"] == "curr_user":
        options["UID"] = pwd.getpwuid(os.getuid()).pw_name
        
    d = getsignal(exp_id, data_name, no_data=no_data, options=options)
    
    return d

def add_coordinate(data_object, new_coordinates, options={}):
    raise NotImplementedError("Coordinate conversions not implemented yet.")