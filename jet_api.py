# -*- coding: utf-8 -*-
"""

.. module:: FLAP/***
   :platform: UNIX
   :synopsis: Signal reading for JET

.. moduleauthor:: M. Vecsei


"""

import os
import pwd
import flap
import numpy as np
from datetime import date

def getsignal(exp_id, source, no_data = False, options={}):
    options_default = {"datapath": flap.config.get("Module JET_API","Datapath")}
    options = {**options_default, **options}
    if options["datapath"].upper() == 'SAL':
        d = getsignal_sal(exp_id, source, no_data=no_data, options=options)
    elif options["datapath"].upper() == 'PPF':
        d = getsignal_ppf(exp_id, source, no_data=no_data, options=options)

    return d

def getsignal_ppf(exp_id, source, no_data = False, options={}):
    """ Signal reading function for the JET API PPF
    Does not require using the FLAP storage
    Input: same as get_data()
    Output: DataObject for the given source
    """
    import ppf
    import getdat as jpf
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
    if os.path.exists(filename):
        return flap.load(filename)

    # obtaining the data from the server
    # if the UID ends with -team, then it looks in the config file, whether
    # there is a team defined under the name given in UID in the config file
    # if so, it will loop over the name of the team members
    split_uid = options["UID"].split("-")
    if split_uid[-1] == "team":
        uid_list = flap.config.get("Module JET_API",options["UID"], evaluate=True)
    else:
        if options["UID"] == "curr_user":
            options["UID"] = pwd.getpwuid(os.getuid()).pw_name
        uid_list = [options["UID"]]

    if split_source[0].upper() == "PPF":
        for uid in uid_list:
            raw_ppf = \
                    ppf.ppfdata(exp_id,split_source[1],split_source[2],
                                seq=options["Sequence"], uid=uid,
                                fix0=options["fix0"], reshape=options["reshape"],
                                no_x=options["no_x"], no_t=options["no_t"],
                                no_data=options["Only Info"])
            [data,x,t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier] = raw_ppf
            if ier == 0:
                break
        if ier != 0:
            raise ValueError("Error reading "+source+" for shot "+str(exp_id)+
                             ", UID "+options["UID"]+":\nIER "+str(ier)+": "+desc)
    elif split_source[0].upper() == "JPF":
        for uid in uid_list:
            raw_jpf = jpf.getdat("/".join(split_source[1:]),
                                 exp_id, nwds = options["nwds"])
            [data, t, nwds, desc, dunits, ier] = raw_jpf
            if ier == 0:
                break
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
    info = "Obtained at "+str(date.today())+", uid "+uid+"\n"
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
        if len(x)>1:
            x_coord = flap.Coordinate(name='X', unit=xunits, values=x,
                                         shape=x.shape,
                                         mode=flap.CoordinateMode(equidistant=False),
                                         dimension_list=[1])
            coordinates.append(x_coord)
    if options["Only Info"] == 0 and (len(data) != 0):
        unit = flap.Unit(name=source,unit=dunits)
        if "x_coord" in locals():
            data = data.reshape((len(t), len(x)))
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
    
def getsignal_sal(exp_id, source, no_data = False, options={}):
    """ Signal reading function for the JET API SAL
    Does not require using the FLAP storage
    Input: same as get_data()
    Output: DataObject for the given source
    """
    # Function used for getting either ppf or jpf data from the server
    #
    from jet.data import sal
    
    options_default = {"Sequence": 0,
                       "UID": "jetppf",
                       "Check Time Equidistant": True,
                       "Cache Data": True}
    options = {**options_default, **options}
    
    #checking if the data is already in the cache
    curr_path = os.path.dirname(os.path.abspath(__file__))
    location = os.path.sep.join(curr_path.split(os.path.sep))
    split_source = source.split("/")
    filename = os.path.join(os.path.sep.join([location,"cached"]),
                            ("_".join(split_source)+"-"+options["UID"]+\
                            "-"+str(exp_id)+".hdf5").lower())
    if os.path.exists(filename):
        return flap.load(filename)

    # obtaining the data from the server
    # if the UID ends with -team, then it looks in the config file, whether
    # there is a team defined under the name given in UID in the config file
    # if so, it will loop over the name of the team members
    split_uid = options["UID"].split("-")
    if split_uid[-1] == "team":
        uid_list = flap.config.get("Module JET_API",options["UID"], evaluate=True)
    else:
        if options["UID"] == "curr_user":
            options["UID"] = pwd.getpwuid(os.getuid()).pw_name
        uid_list = [options["UID"]]

    # has to do some character conversion in the names as SAL uses different
    # naming as the standard ppf names
    character_dict = {'>':'_out_',
                      '<':'_in_',
                      ':':'_sq_',
                      '$':'_do_',
                      ':':'_sq_',
                      '&':'_et_',
                      ' ':'_sp_'}
    node = split_source[2].lower()
    for char in character_dict.keys():
        node = character_dict[char].join(node.split(char))

    signal_error = None
    error_string = str()
    for uid in uid_list:
        try:
            if split_source[0].lower() == "ppf":
                source_string = "/pulse/"+str(exp_id)+"/"+split_source[0].lower()+"/signal/"+uid.lower()+"/"+split_source[1].lower()\
                            +"/"+node
            elif split_source[0].lower() == "jpf":
                source_string = "/pulse/"+str(exp_id)+"/"+split_source[0].lower()+"/"+split_source[1].lower()\
                            +"/"+node+"/data"
            else:
                raise ValueError("Source should be either jpf or ppf")
            if options["Sequence"] != 0:
                source_string = source_string + ":"+str(options["Sequence"])
            raw_signal = sal.get(source_string)
            data = raw_signal.data
            signal_error = raw_signal.error
        except sal.SALException as e:
            error_string = error_string + source_string + " " + str(e) +"\n"
    if not (signal_error is None):
        raise ValueError("Error reading "+source+" for shot "+str(exp_id)+
                         ", UID "+options["UID"]+":\nIER "+str(ier)+": "+desc)
    elif not ("data" in locals()):
        raise ValueError("No data found with errors: \n"+error_string)
    # converting data to flap DataObject
    # possible names for the time coordinate
    time_names=['time', 't', 'jpf time vector', 'ppf time vector', 'time vector', 'the ppf t-vector.']
    coordinates = []
    data = data.reshape(raw_signal.shape)
    info = "Obtained at "+str(date.today())+", uid "+uid+"\n"
    coord_dimension = 0
    for coord in raw_signal.dimensions:
        name =  coord.description
        values = np.array(coord.data)
        unit = coord.units
        equidist = False
        if name.lower() in time_names or coord.temporal is True:
            tunit = ["secs", "s", "seconds", "second"]
            if unit.lower() in tunit:
                unit = "Second"
            name = 'Time'
            equidist = False
            if options['Check Time Equidistant'] is True and (unit == "Second" or coord.temporal is True)\
               and (len(values)==1 and values[0]==-1):
                timesteps = np.array(values[1:])-np.array(values[0:-1])
                equidistance = np.linalg.norm(timesteps-timesteps[0])/np.linalg.norm(timesteps)
                if equidistance < 1e-6:
                    info = info + "Time variable is taken as equidistant to an accuracy of "+\
                           str(equidistance)+"\n"
                    coord_object = flap.Coordinate(name=name, unit=values, start = values[0],
                                                   shape=values.shape, step=np.mean(timesteps),
                                                   mode=flap.CoordinateMode(equidistant=True),
                                                   dimension_list=[0])
                    equidist = True
        if equidist is False:
            coord_object = flap.Coordinate(name=name, unit=unit, values=values,
                                           shape=values.shape,
                                           mode=flap.CoordinateMode(equidistant=False),
                                           dimension_list=[coord_dimension])
        coordinates.append(coord_object)
        coord_dimension = coord_dimension+1

    unit = flap.Unit(name=source,unit=raw_signal.units)
    signal = flap.DataObject(data_array=data, data_unit=unit,
                             coordinates=coordinates, exp_id=str(exp_id),
                             data_title=raw_signal.description, data_shape=data.shape,
                             info=info)
    if "summary" in dir(raw_signal):
        signal.info = signal.info+raw_signal.summary().description+"\n"
    if options["Cache Data"] is True:
        flap.save(signal, filename)
    return signal


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
               UID can also beam a "-team" in this case it has to be defined in the flap_defaults.cfg
               under the JET_API Module
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
            "nwds": 0,
            "datapath": flap.config.get("Module JET_API","Datapath")}
    options = {**options_default, **options}
    
    d = getsignal(exp_id, data_name, no_data=no_data, options=options)

    return d

def add_coordinate(data_object, new_coordinates, options={}):
    raise NotImplementedError("Coordinate conversions not implemented yet.")