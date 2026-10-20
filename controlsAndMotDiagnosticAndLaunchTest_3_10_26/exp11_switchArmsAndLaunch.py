# Experiment: MOT + MRR Launch with PGC
#
# Purpose: Form MOT on MOT arm, accelerate MRR stage, PGC cool,
# switch to MRR arm, open shutter, PGC in moving frame, image.
#
# Globals:
#   t_load:             MOT loading time in seconds (default 4.0)
#   t_dark:             shutter-closed dark time in seconds (default 0.028)
#   t_hold:             time between molasses images in seconds (default 0.02)
#   t_redshift:         PGC duration before closing shutter in seconds (default 0.005)
#   t_redshift_launch:     PGC duration in moving frame in seconds (default 0.005)
#   v_laser_jump_rel:   red shift amplitude in volts (default -0.05)
#   v_stage:            stage velocity in mm/s (default 5)
#   a_stage:            stage acceleration in mm/s^2 (default 10000)
#   xf_stage:           stage move distance in mm (default 1)

import numpy as np
from labscript import *
from labscriptlib.eedmLabscriptOps.connection_table_eedmW52 import ct
from labscript import compiler


def set_beam_path(t, mode):
    """Switch beam path between 'MOT' and 'MRR'."""
    if mode == 'MOT': #set LCR to V2 (low) #note, old liquid crystal controller is BACKWARDS from new one, high = MOT 
        #MOT_SHUTTER_do.open(t)
        #MRR_SHUTTER_do.close(t)
        LCR_do.go_high(t)
        LCR_SPLIT_do.go_low(t)
        #LCR_HOR_ao.constant(t, 0.0)
        #LCR_TOP_ao.constant(t, 0.0)
    elif mode == 'MRR': #set LCR to V1 (high)
        #MOT_SHUTTER_do.close(t)
        #MRR_SHUTTER_do.open(t)
        LCR_do.go_low(t)
        LCR_SPLIT_do.go_high(t)
        #LCR_HOR_ao.constant(t, 5.0)
        #LCR_TOP_ao.constant(t, 5.0)
    else:
        raise ValueError(f"Unknown beam path mode: {mode}. Use 'MOT' or 'MRR'.")


if __name__ == "__main__":
    ct()
    start()

    t = 0

    # --- Runmanager globals with defaults ---
    try:
        t_load
    except NameError:
        t_load = 4.0

    try:
        t_dark
    except NameError:
        t_dark = 0.028

    try:
        t_hold
    except NameError:
        t_hold = 0.02

    try:
        t_redshift
    except NameError:
        t_redshift = 0.005

    try:
        t_redshift_launch
    except NameError:
        t_redshift_launch = 0.005

    try:
        v_laser_jump_rel
    except NameError:
        v_laser_jump_rel = -0.05

    try:
        v_stage
    except NameError:
        v_stage = 5

    try:
        a_stage
    except NameError:
        a_stage = 10000

    try:
        xf_stage
    except NameError:
        xf_stage = 1

    # --- Stage kinematics ---
    MRR_shutter_min_time = 0.012
    MOT_coil_delay = 0.003

    t_accel = v_stage / a_stage
    d_accel = 0.5 * a_stage * t_accel**2
    xi_stage = np.round((150 - d_accel - 0.25), 1)

    if d_accel > abs(xf_stage - xi_stage):
        raise ValueError("Target velocity too high for move distance and acceleration.")

    # --- Configure stage ---
    bbd301.set_velocity(1, v_stage)
    bbd301.set_acceleration(1, a_stage)
    bbd301.set_move_distance(1, xf_stage)
    bbd301.set_reset_position(1, MRR_START_POS)
    bbd301.set_start_position(1, MRR_START_POS)

    if 0: #print stage config for debugging
        print(f"Stage config:")
        print(f"  xi_stage:    {xi_stage} mm")
        print(f"  xf_stage:    {xf_stage} mm")
        print(f"  v_stage:     {v_stage} mm/s")
        print(f"  a_stage:     {a_stage} mm/s^2")
        print(f"  t_accel:     {t_accel*1000:.3f} ms")
        print(f"  d_accel:     {d_accel:.4f} mm")
        print(f"  stage_delay: {MRR_TRIG_DELAY*1000:.1f} ms")
        print(f"  t_stage_ready: {(t_accel + MRR_TRIG_DELAY)*1000:.3f} ms")
        t_move_total = xf_stage / v_stage  # approximate total move time at constant v
        print(f"  t_move_total (approx): {t_move_total*1000:.1f} ms")
        t_experiment = t_load + 0.01 + MOT_coil_delay + t_redshift + max(t_dark, t_accel + MRR_TRIG_DELAY) + 0.005 + 0.001 + 50e-6 + t_redshift_launch + 0.001 + 4*t_hold + 0.01 + 2
        print(f"  t_experiment (approx): {t_experiment*1000:.1f} ms")
        if t_move_total > t_experiment:
            print(f"  WARNING: stage may still be moving at end of experiment!")
        else:
            print(f"  OK: stage should finish before experiment ends.")
    if 1: #Stage position at shutter open, TTL->motion jitter 0.1-6 ms
        def stage_pos(dt):
            if dt <= 0:
                return 0.0, 0.0
            if dt < t_accel:
                return 0.5*a_stage*dt**2, a_stage*dt
            return d_accel + v_stage*(dt - t_accel)+MRR_START_POS, v_stage

        print("Position at MRR shutter open (t_shutter = TTL + MRR_TRIG_DELAY):")
        for jitter in (0.0001, 0.006):
            dt = MRR_TRIG_DELAY - jitter
            x, v = stage_pos(dt)
            print(f"  jitter {jitter*1e3:5.1f} ms -> x = {x:.4f} mm")

    # ============================================================
    # Initialize outputs -- MOT arm, everything normal
    # ============================================================
    MOT_COIL_do.go_low(t)
    REPUMP_SHUTTER_do.open(t)
    MAIN_REL_JUMP_do.go_low(t)
    #REPUMP_REL_JUMP_do.go_low(t)
    MRR_TRIG_do.go_low(t)
    #MAIN_JUMP_AMP_ao.constant(t, v_laser_jump_rel)
    MAIN_JUMP_AMP_ao.constant(t, 0.0)
    set_beam_path(t, 'MOT')
    MOT_SHUTTER_do.close(t)
    MRR_SHUTTER_2_do.close(t)
    MRR_SHUTTER_do.open(t)
    PROBE_SHUTTER_do.close(t)
    

    #LCR_HOR_ao.constant(t, 5.0)
    #LCR_TOP_ao.constant(t,0.886)
    
    #MOT_SHUTTER_do.open(t)
    #MRR_SHUTTER_do.open(t)
    #LCR_do.go_high(t)
    #LCR_HOR_ao.constant(t, 5.0)
    #LCR_BOT_do.go_low(t)
    #LCR_TOP_ao.constant(t, 5.0)

    t += 0.25
    #MOT_SHUTTER_do.open(t)
    MOT_COIL_do.go_high(t)
    MOT_SHUTTER_do.open(t)


    # ============================================================
    # Load MOT on MOT arm
    # ============================================================
    add_time_marker(t, "Loading MOT", verbose=True)
    t += t_load

    # ============================================================
    # Reference image -- MOT in home position
    # ============================================================
    #add_time_marker(t, "Reference image", verbose=True)
    #my_ids_camera.expose(
    #    t=t, name='reference', frametype='atom',
    #    trigger_duration=1*ms
    #)
    t += 0.01

    # ============================================================
    # MOT coils off for molasses
    # ============================================================
    add_time_marker(t, "MOT coils off", verbose=True)
    MOT_COIL_do.go_low(t)
    t += 0.01

    # ============================================================
    # PGC on MOT arm -- red shift on
    # ============================================================
    add_time_marker(t, "PGC on MOT arm", verbose=True)
    MAIN_JUMP_AMP_ao.constant(t-0.001, v_laser_jump_rel)
    MAIN_REL_JUMP_do.go_high(t)
    t += t_redshift

    # ============================================================
    # Close MOT shutter -- dark time begins
    # Atoms in the dark, stage still accelerating
    # ============================================================
    add_time_marker(t, "Close MOT shutter -- dark time", verbose=True)
    MOT_SHUTTER_do.close(t)
    
    #scope_trig_do.go_high(t)#using this as a scope trigger. Delete this
    #scope_trig_do.go_low(t+0.1)
    #t_trig = t

    t += 0.002
    #REPUMP_SHUTTER_do.close(t)
    MAIN_REL_JUMP_do.go_low(t)
    MAIN_JUMP_AMP_ao.constant(t+0.001, 0.00)
    

    add_time_marker(t, "Switch to MRR arm", verbose=True)
    set_beam_path(t, 'MRR')
    #LCR_TOP_ao.constant(t, 5.0)

    t += t_dark 
    # ============================================================
    # Switch to MRR arm -- LCRs flip during dark time
    # ============================================================

    MRR_SHUTTER_2_do.open(t)
        
    scope_trig_do.go_high(t)#using this as a scope trigger. Delete this
    scope_trig_do.go_low(t+0.3)

    # ============================================================
    # Trigger stage to start accelerating
    # Do this early so it reaches speed during PGC + dark time
    # ============================================================
    print('v stage =', v_stage)
    print('a stage =', a_stage)
    print('v/a =', v_stage/a_stage)
    add_time_marker(t-(MRR_TRIG_DELAY), "Trigger MRR stage", verbose=True)
    if MRR_TRIG_BOOL:
        MRR_TRIG_do.go_high(t-(MRR_TRIG_DELAY)) #trigger at time such that stage starts moving at t=0
        MRR_TRIG_do.go_low(t-(MRR_TRIG_DELAY) + 0.011)
        #scope_trig_do.go_high(t-(v_stage/a_stage + 0.05))#using this as a scope trigger. Delete this
        #scope_trig_do.go_low(t-(v_stage/a_stage + 0.05) + 0.011)
        #ai1.acquire(label='MRR_position', start_time = t-(v_stage/a_stage + 0.05), end_time =t-(v_stage/a_stage + 0.05)+0.2)

    # ============================================================
    # Wait for photodiode trigger confirming shutter is open
    # ============================================================
    add_time_marker(t, "Wait for shutter open trigger", verbose=True)
    wait("mrr_shutter_open", t, timeout=0.1)
    
    

    
    #============================================================
    #PGC in moving frame -- red shift on ~50 us after shutter open
    #============================================================
    
    #MAIN_JUMP_AMP_ao.constant(t, -0.1)

    #my_ids_camera.expose(
    #t=t+0.0005, name=f'molasses', frametype='atom',
    #trigger_duration=1*ms)

    
    #t += 1000e-6 #works for drop
    t += 0.0035


    add_time_marker(t, "PGC in moving frame", verbose=True)
    MAIN_REL_JUMP_do.go_high(t)
    MAIN_JUMP_AMP_ao.ramp(t-0.000, 0.005, -0.0, -.15, 1e5) #(t, duration, initial, final, samplerate)
    #t += 0.01 #works for drop
    t+= 0.0035
    MRR_SHUTTER_do.close(t)
    REPUMP_SHUTTER_do.close(t+0.01) 
    drop_time = t  # atoms begin free fall here
    compiler.shot_properties['drop_time'] = drop_time
    photon_counter.acquire(t=drop_time+t_PMT_delay, number_of_counts=int(.15e5)) #bins are 10us, hard coded
    ai0.acquire(label='TOF_florescence', start_time = drop_time+t_PMT_delay, end_time =drop_time+t_PMT_delay+0.2)
    ai1.acquire(label='RelJumpModulation', start_time = drop_time+t_PMT_delay, end_time = drop_time+t_PMT_delay+0.2)
    
    #scope_trig_do.go_high(drop_time+t_PMT_delay)
    #scope_trig_do.go_low(drop_time+t_PMT_delay+0.02)
    t_trig = t



    #============================================================
     #Jump back on resonance for imaging
    #============================================================
    t += 0.002
    MAIN_REL_JUMP_do.go_low(t)
    
    #jump probe to blue
    t += 0.001
    #MAIN_JUMP_AMP_ao.constant(t, 0.006)

    t += t_PMT_delay-0.01
    
    #REPUMP_SHUTTER_do.open(t)
    PROBE_SHUTTER_do.open(t)

    t += 0.05

    MAIN_JUMP_AMP_ao.sine(t, 0.1, 0.01 , 2*np.pi*6.123e3,0,0.006,1e5) #time, duration, amp, angfreq, phase, dc_offset samplerate
    add_time_marker(t, "Jump to blue", verbose=True)
    MAIN_REL_JUMP_do.go_high(t)
    MAIN_REL_JUMP_do.go_low(t+0.07)
    '''
    t += t_camera_delay
    my_ids_camera.expose(
    t=t, name=f'molasses', frametype='atom',
    trigger_duration=1*ms)'''
    '''
    # ============================================================
    # Molasses images
    # ============================================================
    for i in range(1):
        add_time_marker(t, f"Molasses image {i}", verbose=True)
        my_ids_camera.expose(
            t=t, name=f'molasses_{i}', frametype='atom',
            trigger_duration=1*ms
        )
        t += t_hold
    '''
    # ============================================================
    # Cleanup -- return to MOT arm, normal state
    # ============================================================
    t += .2
    MAIN_REL_JUMP_do.go_low(t)
    #REPUMP_REL_JUMP_do.go_low(t)
    MOT_COIL_do.go_high(t)
    set_beam_path(t, 'MOT')
    REPUMP_SHUTTER_do.open(t)
    MOT_SHUTTER_do.open(t)
    MRR_SHUTTER_2_do.close(t)

    t += 2
    stop(t)