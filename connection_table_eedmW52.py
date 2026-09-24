
from labscript import *
from labscript import AnalogOut
from labscript_devices.DummyPseudoclock.labscript_devices import DummyPseudoclock
from labscript_devices.DummyIntermediateDevice import DummyIntermediateDevice
from labscript_devices.NI_DAQmx.models.NI_PXIe_6738 import NI_PXIe_6738
from labscript_devices.NI_DAQmx.models.NI_PXIe_6363 import NI_PXIe_6363
from labscript_devices.IMAQdxCamera.labscript_devices import IMAQdxCamera
from labscript_devices.PrawnBlaster.labscript_devices import PrawnBlaster
from user_devices.NI_DAQmx_PhotonCounter.labscript_devices import NI_DAQmxPhotonCounter

#from labscript_devices import IDS_PeakCamera
from user_devices.IDS_PeakCamera.labscript_devices import IDS_PeakCamera
from user_devices.BBD301.labscript_devices import BBD301

from labscript.remote import RemoteBLACS

def ct():
    # Remote worker host: everything below runs its BLACS worker on ni-eedm private ip
    pxi = RemoteBLACS(name='pxi_crate', host='eedm-ni.local', port=7341)

    PrawnBlaster(
        name                = 'prawn', 
        com_port            = 'COM11', 
        num_pseudoclocks    = 1,
        pico_board          = 'pico2'
    )

    prawn.set_property("worker", pxi.name, "connection_table_properties") #consider moving to w52

    SLOT_6738 = "PXI1Slot3"   # EXACT name from NI-MAX
    SLOT_6363 = "PXI1Slot2"   # EXACT name from NI-MAX

    NI_PXIe_6738(
        name                    = 'NI6738',
        parent_device           = prawn.clocklines[0],
        MAX_name                = SLOT_6738,
        clock_terminal          = '/PXI1Slot3/PFI0',         # <-- empty string means "use internal/Onboard clock"
        max_AO_sample_rate      = 100000.00000000001,
        worker=pxi
    )

    NI_PXIe_6363(
        name                    = 'NI6363',
        parent_device           = prawn.clocklines[0],
        MAX_name=SLOT_6363,
        clock_terminal          = '/PXI1Slot2/PFI12',         # <-- empty string means "use internal/Onboard clock"
        acquisition_rate=100e3,
        worker=pxi
    )

    # Define digital outs on 6363
    #make sure to define an even number of channels for Labscript to be happy
    Shutter(name='REPUMP_SHUTTER_do', parent_device=NI6363, connection='port0/line0',
        delay=(5.94e-3, 5.6e-3), open_state=1) #Delays coppied from MRR, not directly tested!
    DigitalOut(name='MAIN_REL_JUMP_do', parent_device=NI6363, connection='port0/line1')
    DigitalOut(name='LCR_SPLIT_do', parent_device=NI6363, connection='port0/line11')
    DigitalOut(name='MOT_COIL_do', parent_device=NI6363, connection='port0/line3')
    DigitalOut(name='MRR_TRIG_do', parent_device=NI6363, connection='port0/line4')
    Shutter(name='MOT_SHUTTER_do', parent_device=NI6363, connection='port0/line5',
            delay=(6e-3, 6e-3), open_state=1) #delays are a guess.
    Shutter(name='MRR_SHUTTER_do', parent_device=NI6363, connection='port0/line6',
            delay=(5.94e-3, 5.6e-3), open_state=1) #Delays measured and set 2026-FEB-26. See lab notebook for details.
    DigitalOut(name='LCR_do', parent_device=NI6363, connection='port0/line7')
    DigitalOut(name='scope_trig_do', parent_device=NI6363, connection='port0/line8')
    Shutter(name='MRR_SHUTTER_2_do', parent_device=NI6363, connection='port0/line9',
            delay=(5.47e-3, 7.08e-3), open_state=1) #Delays measured and set 2026-MAY-04. See lab notebook for details.
    Shutter(name='PROBE_SHUTTER_do', parent_device=NI6363, connection='port0/line10',
            delay=(6e-3, 6e-3), open_state=1) #Delays are a guess

    DigitalOut(name='dummy_do', parent_device=NI6363, connection='port0/line28')    

    # Define analog outs on 6363. Using these as digital outs for now.
    AnalogOut(name='LCR_HOR_ao', parent_device=NI6363, connection='ao0')
    AnalogOut(name='LCR_TOP_ao', parent_device=NI6363, connection='ao1')

    # Define analog outs on 6738
    # Make sure to define an even number of channels for Labscript to be happy
    AnalogOut(name='ao0', parent_device=NI6738, connection='ao0')
    AnalogOut(name='ao1', parent_device=NI6738, connection='ao1')
    AnalogOut(name='ao2', parent_device=NI6738, connection='ao2')
    AnalogOut(name='MAIN_JUMP_AMP_ao', parent_device=NI6738, connection='ao3')
    AnalogOut(name='REPUMP_JUMP_AMP_ao', parent_device=NI6738, connection='ao4')
    AnalogOut(name='ao5', parent_device=NI6738, connection='ao5')
    AnalogOut(name='ao6', parent_device=NI6738, connection='ao6')
    AnalogOut(name='ao7', parent_device=NI6738, connection='ao7')

    # Define analog ins on 6363
    AnalogIn(name='ai0', parent_device=NI6363, connection='ai0')
    AnalogIn(name='ai1', parent_device=NI6363, connection='ai1')

    # Again, make sure to define an even number of channels for Labscript to be happy
    #do0 = DigitalOut('do0', ao6738, 'port0/line1')
    DigitalOut('do1_dummy', NI6738, 'port0/line0')  # dummy to make even number of digital lines

    # cam_serial_number2 = '4102766024' #sn yellow 4103389953 blue 4108596607 red 4108596608 green 4102766024
    cam = IDS_PeakCamera( #blue
        name='my_ids_camera',
        parent_device=NI6738,
        connection='port0/line1', 
        serial_number='4108596607', #sn yellow 4103389953 blue 4108596607 red 4108596608
        minimum_recovery_time=0.005, # Override the default
        trigger_edge_type='rising', # We will trigger via software
        trigger_duration=0.05,
        camera_attributes={
            'trigger': 'On', # On/Off
            'format': 'Mono12', # Mono8/Mono12
            'exposure': 1, # global property exposure_time overrides this, but we need to set it here too for the BLACS preview mode
            'fps': 5.0,   # required by base class (will be skipped at runtime)
            'gain': 0.0 #gloabal gain overides this, but we need to set it here too for the BLACS preview mode,
        },
        manual_mode_camera_attributes={ # BLACS preview mode
            "trigger": "Off",
            "format": "Mono8",
            "exposure": 2.5,
            "fps": 5.0,               # shown in BLACS, not applied to hardware
            "gain": 1.0,
        },
    )


    '''
    cam = IDS_PeakCamera( #green
         name=f'my_ids_camera',
         parent_device=NI6738, 
         connection='port0/line1', #cheat: fake connection to make labscript happy
         serial_number=4102766024, 
         minimum_recovery_time=0.005, # Override the default
         trigger_edge_type='rising', 
         trigger_duration=100*ms,
         camera_attributes={
             'trigger': 'On', # On/Off
             'format': 'Mono8', # Mono8/Mono12
             'exposure': 9.0, # 9 ms, using the wrapper property
             'fps': 5.0,   # required by base class (will be skipped at runtime)
             'gain': 1.0,
             'roi': [0,0,1440,1080],
             'binning': [1,1]
         },
         manual_mode_camera_attributes={ # BLACS preview mode
             "trigger": "Off",
             "format": "Mono8",
             "exposure": 20.0,
             "fps": 5.0,               # shown in BLACS, not applied to hardware
             "gain": 1.0,
             "roi": [0,0,1440,1080],
             "binning": [1,1]
         },
     )'''

    bbd = BBD301(
        name='bbd301',
        parent_device=prawn.clocklines[0],
        serial_number='103512594',  # your actual serial number
        num_channels=1,
        worker=pxi
    )

    photon_counter = NI_DAQmxPhotonCounter(
        name='photon_counter',
        parent_device=prawn.clocklines[0],
        #connection='internal',
        MAX_name=SLOT_6363,
        counter_channel='ctr3',
        photon_input_terminal=f'/{SLOT_6363}/PFI5',  # wire your photon source here
        sample_clock_terminal='',  # None = internal clock (matches DummyPseudoclock)
        counter_sample_rate=100000,
        start_trigger_terminal=f'/{SLOT_6363}/PFI12',
        worker=pxi
    )



if __name__ == '__main__':
    ct()
    start()
    stop(0.001)