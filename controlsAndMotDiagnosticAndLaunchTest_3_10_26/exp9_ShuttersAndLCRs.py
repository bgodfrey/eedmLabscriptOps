# Manual beam line test
# 
# Purpose: Static test of all 6 TTL lines. Set each flag below to 1 or 0
# to control which outputs go high when the coil trigger fires.
#
# The coil trigger (MOT_COIL_do) goes high at t = 10 ms.
# All other outputs are set at the same time as the coil trigger.

from labscript import *
from labscriptlib.eedmLabscriptOps.connection_table_eedmW52 import ct

if __name__ == "__main__":
    ct()
    start()

    # ============================================================
    # Set these to 1 or 0 to enable/disable each output. Connection table is pretty crazy
    # ============================================================
    REPUMP_SHUTTER    = 1   # port0/line0 - repump shutter
    MAIN_REL_JUMP     = 0   # port0/line1 - main rel jump
    REPUMP_REL_JUMP   = 0   # port0/line2 - repump rel jump
    MRR_TRIG          = 0   # port0/line4 - MRR trigger
    MOT_SHUTTER       = 1   # port0/line5 - MOT shutter
    MRR_SHUTTER       = 0   # port0/line6 - MRR shutter
    LCR               = 0   # port0/line7 - LCR waveplate
    # ============================================================

    t = 0

    # Everything starts low / closed
    MOT_COIL_do.go_low(t)
    REPUMP_SHUTTER_do.close(t)
    MAIN_REL_JUMP_do.go_low(t)
    #REPUMP_REL_JUMP_do.go_low(t)
    MRR_TRIG_do.go_low(t)
    MOT_SHUTTER_do.close(t)
    MRR_SHUTTER_do.close(t)
    LCR_do.go_low(t)

    # Coil trigger fires, outputs go to configured state
    t = 10e-3

    MOT_COIL_do.go_high(t)

    if REPUMP_SHUTTER:
        REPUMP_SHUTTER_do.open(t)
    if MAIN_REL_JUMP:
        MAIN_REL_JUMP_do.go_high(t)
    #if REPUMP_REL_JUMP:
    #    REPUMP_REL_JUMP_do.go_high(t)
    if MRR_TRIG:
        MRR_TRIG_do.go_high(t)
    if MOT_SHUTTER:
        MOT_SHUTTER_do.open(t)
    if MRR_SHUTTER:
        MRR_SHUTTER_do.open(t)
    if LCR:
        LCR_do.go_high(t)

    t += 0.1
    stop(t)