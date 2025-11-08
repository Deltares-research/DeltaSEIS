"""
Top Mute Method - Complete Documentation and Examples
=====================================================

The top_mute() method zeros out seismic data above a specified horizon,
commonly used to remove water column noise or data above first breaks.

PARAMETERS:
-----------
horizon : float or array
    - float: constant mute time (ms) for all traces
    - array: time value (ms) per trace (must match number of traces)

shift_ms : float (default=0.0)
    Shift the mute boundary in milliseconds
    - Positive: mute MORE data (extend mute zone DEEPER)
    - Negative: mute LESS data (preserve MORE signal ABOVE horizon)
    Example: shift_ms=-1.0 means mute everything above (horizon - 1 ms)

taper_ms : float (default=0.0)
    Cosine taper length in milliseconds at mute boundary
    - 0: hard cutoff (instant transition from 0 to 1)
    - >0: smooth ramp from 0 (muted) to 1 (preserved)
    Recommended: 2-5 ms to avoid edge artifacts

inplace : bool (default=True)
    - True: modifies self.data in place
    - False: returns muted data without modifying original

COMMON USE CASES:
-----------------
"""

import numpy as np
import deltaseis
from pathlib import Path

# Example data loading
segy_folder = Path(r'D:\Projects\DIS_Offshore\data\xstar')
segy_paths = list(segy_folder.glob('*.sgy'))
if len(segy_paths) > 0:
    s = deltaseis.Segy_edit(segy_paths[0])
    s.xstar_split('full')
    s.get_seabed_pick(10, 100, 9, 3, truncate=10)
    s.filter_horizon_savgol('seabed_pick', 'seabed_pick_savgol', 501, 4)
    
    seis = deltaseis.Seismic(np.array(s.trace_data).T, fs=50_000, dx=0.4)
    seis.time_power_gain(3)

print("="*70)
print("EXAMPLE 1: Basic Usage - Mute at Seabed")
print("="*70)
print("""
# Mute exactly at seabed pick (hard cutoff)
seis.top_mute(s.seabed_pick, shift_ms=0.0, taper_ms=0.0)

USE CASE: Remove all data above seabed
CAUTION: Hard cutoff may cause ringing artifacts
""")

print("="*70)
print("EXAMPLE 2: Conservative Mute - Preserve Data Above Seabed")
print("="*70)
print("""
# Mute 1-2 ms ABOVE seabed (preserve near-seabed data for QC)
seis.top_mute(s.seabed_pick, shift_ms=-1.5, taper_ms=0.0)

USE CASE: Keep some data above seabed for quality control
BEST PRACTICE: Use negative shift to be conservative
TIP: Check if seabed pick is accurate before aggressive muting
""")

print("="*70)
print("EXAMPLE 3: Smooth Mute with Taper (RECOMMENDED)")
print("="*70)
print("""
# Mute 1 ms above seabed with 2 ms cosine taper
seis.top_mute(s.seabed_pick, shift_ms=-1.0, taper_ms=2.0)

USE CASE: Production processing - balances noise removal with artifact prevention
WHY TAPER: Avoids sharp discontinuities that can cause ringing
TAPER LENGTH: 2-5 ms is typical (depends on dominant frequency)
""")

print("="*70)
print("EXAMPLE 4: Aggressive Mute - Remove Noisy Near-Seabed Data")
print("="*70)
print("""
# Mute 2 ms BELOW seabed (extend mute zone deeper)
seis.top_mute(s.seabed_pick, shift_ms=+2.0, taper_ms=3.0)

USE CASE: Remove ringing or multiples near seabed
CAUTION: You're removing signal data - use sparingly
TIP: Only use positive shift if confident about noise zone extent
""")

print("="*70)
print("EXAMPLE 5: Using Filtered Horizon")
print("="*70)
print("""
# Use smoothed seabed pick to avoid mute time fluctuations
seis.top_mute(s.seabed_pick_savgol, shift_ms=-1.0, taper_ms=2.0)

USE CASE: When seabed pick has high-frequency noise
ADVANTAGE: Smoother mute boundary, less trace-to-trace variation
BEST PRACTICE: Apply Savitzky-Golay filter to horizon before muting
""")

print("="*70)
print("EXAMPLE 6: Constant Time Mute")
print("="*70)
print("""
# Mute above a fixed time (e.g., 25 ms for all traces)
seis.top_mute(25.0, shift_ms=0.0, taper_ms=3.0)

USE CASE: Land data with constant first break time
ALSO: Initial QC to remove very early noise
TIP: Float value applies same time to all traces
""")

print("="*70)
print("EXAMPLE 7: Workflow Integration")
print("="*70)
print("""
# Typical processing workflow:
s = deltaseis.Segy_edit('data.sgy')
s.xstar_split('full')
s.get_seabed_pick(10, 100, 9, 3)
s.filter_horizon_savgol('seabed_pick', 'seabed_pick_savgol', 501, 4)

# Transfer to Seismic object
seis = deltaseis.Seismic(np.array(s.trace_data).T, fs=50_000, dx=0.4)
seis.time_power_gain(3)

# Apply deconvolution (amplifies noise!)
seis.signature_deconvolution(6414, 30.5, 31.7, method='wiener', epsilon=0.05)

# NOW apply top mute to remove amplified noise above seabed
seis.top_mute(s.seabed_pick_savgol, shift_ms=-1.0, taper_ms=2.0)

# Bandpass filter to further reduce noise
seis.bandpass_filter(lowcut=1000, highcut=20000)

# Write back
s.trace_data = seis.data.T
s.write('output.sgy')
""")

print("="*70)
print("EXAMPLE 8: Return Modified Data (Non-destructive)")
print("="*70)
print("""
# Get muted data without modifying original
muted_data = seis.top_mute(s.seabed_pick, shift_ms=-1.0, 
                           taper_ms=2.0, inplace=False)

# Now you can compare
difference = seis.data - muted_data  # Shows what was removed

USE CASE: QC and comparison
ADVANTAGE: Keep original data intact for testing
""")

print("="*70)
print("PARAMETER SELECTION GUIDE")
print("="*70)
print("""
shift_ms SELECTION:
------------------
  -2.0 ms: Very conservative, maximum signal preservation
  -1.0 ms: Recommended for most cases
   0.0 ms: Mute exactly at horizon
  +1.0 ms: Remove potentially noisy data near horizon
  +2.0 ms: Aggressive mute (use with caution)

taper_ms SELECTION:
------------------
   0.0 ms: No taper (hard cutoff) - use only if necessary
   1.0 ms: Short taper for high-frequency data (>10 kHz)
   2.0 ms: Standard taper for typical seismic (RECOMMENDED)
   3.0 ms: Long taper for low-frequency data (<5 kHz)
   5.0 ms: Very long taper to minimize artifacts

RULE OF THUMB: taper_ms ≈ 1-2 periods of dominant frequency
Example: For 2 kHz data (0.5 ms period), use 1-2 ms taper
""")

print("="*70)
print("TROUBLESHOOTING")
print("="*70)
print("""
PROBLEM: Seeing ringing artifacts after mute
SOLUTION: Increase taper_ms (try 3-5 ms)

PROBLEM: Losing too much signal data
SOLUTION: Use more negative shift_ms (e.g., -2.0 instead of -1.0)

PROBLEM: Still seeing noise above seabed after mute
SOLUTION: 
  1. Check if horizon is accurate (plot it!)
  2. Use less negative or zero shift_ms
  3. Consider using smoothed horizon (seabed_pick_savgol)

PROBLEM: ValueError about horizon length
SOLUTION: Ensure horizon array length matches number of traces in seis.data

PROBLEM: NaN values in horizon causing issues
SOLUTION: Method handles NaNs automatically - those traces won't be muted
          Check warning message for NaN count
""")

print("\n" + "="*70)
print("COMPLETE WORKING EXAMPLE")
print("="*70)

if len(segy_paths) > 0:
    print("\nRunning complete example...")
    
    # Fresh copy for demonstration
    seis_demo = deltaseis.Seismic(np.array(s.trace_data).T[:, 5000:5100], fs=50_000, dx=0.4)
    seis_demo.time_power_gain(3)
    
    print("\nOriginal data shape:", seis_demo.data.shape)
    print("Horizon shape:", s.seabed_pick_savgol[5000:5100].shape)
    
    # Apply mute with recommended parameters
    seis_demo.top_mute(s.seabed_pick_savgol[5000:5100], 
                       shift_ms=-1.0, 
                       taper_ms=2.0)
    
    print("\n✓ Top mute successfully applied!")
    print("\nRecommended next steps:")
    print("  1. Visual QC: Check seismic section to verify mute quality")
    print("  2. Horizon QC: Verify seabed pick accuracy before muting")
    print("  3. Parameter tuning: Adjust shift_ms and taper_ms as needed")
    print("  4. Integration: Add to your processing workflow")

print("\n" + "="*70)
print("For more information, see: deltaseis.base_seismic.Seismic.top_mute")
print("="*70)
