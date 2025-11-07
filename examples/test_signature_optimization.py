"""
Test automatic signature window optimization
"""
import numpy as np
import deltaseis
from pathlib import Path

# Load data
segy_folder = Path(r'D:\Projects\DIS_Offshore\data\xstar')
segy_paths = [f for f in segy_folder.iterdir() if f.suffix.lower() in {'.sgy', '.seg', '.segy'}]
segy_path = segy_paths[0]

print(f"Loading: {segy_path.name}\n")
segy = deltaseis.Segy_edit(segy_path)
segy.xstar_split('full')

trace_array = np.array(segy.trace_data)
if trace_array.shape[0] == 12404:
    trace_array = trace_array.T

# Test on a subset for speed
trace_subset = trace_array[:, 6000:6800]

fs = 50_000
seis = deltaseis.Seismic(trace_subset.copy(), fs=fs, dx=0.4)
seis.time_power_gain(3)
seis.trace_averaging(1)

# Original parameters (user-specified)
trace_number = 415  # 6415 - 6000 (subset offset)
start_time_ms = 30.5
end_time_ms = 31.7

print("="*60)
print("TESTING AUTOMATIC SIGNATURE OPTIMIZATION")
print("="*60)

# Test different optimization criteria
criteria = ['spectral_flatness', 'energy', 'kurtosis', 'bandwidth']

results = {}
for criterion in criteria:
    print(f"\n{'='*60}")
    print(f"Testing criterion: {criterion.upper()}")
    print(f"{'='*60}")
    
    result = seis.optimize_signature_window(
        trace_number=trace_number,
        start_time_ms=start_time_ms,
        end_time_ms=end_time_ms,
        time_search_ms=2.0,  # Search ±2 ms around center
        length_search_ms=0.5,  # Adjust length by ±0.5 ms
        criterion=criterion,
        plot=True  # Generate diagnostic plots
    )
    
    results[criterion] = result
    
    # Show detailed results
    print(f"\n  Original window:")
    print(f"    Start: {start_time_ms:.3f} ms")
    print(f"    End:   {end_time_ms:.3f} ms")
    print(f"    Length: {end_time_ms - start_time_ms:.3f} ms")
    
    print(f"\n  Optimized window:")
    print(f"    Start: {result['start_time_ms']:.3f} ms")
    print(f"    End:   {result['end_time_ms']:.3f} ms")
    print(f"    Length: {result['end_time_ms'] - result['start_time_ms']:.3f} ms")
    
    print(f"\n  Changes:")
    print(f"    Start shift: {result['start_time_ms'] - start_time_ms:+.3f} ms")
    print(f"    End shift:   {result['end_time_ms'] - end_time_ms:+.3f} ms")
    print(f"    Quality score: {result['quality_score']:.6f}")

print("\n" + "="*60)
print("SUMMARY OF ALL CRITERIA")
print("="*60)

for criterion, result in results.items():
    print(f"\n{criterion.upper()}:")
    print(f"  Window: {result['start_time_ms']:.2f} - {result['end_time_ms']:.2f} ms")
    print(f"  Score:  {result['quality_score']:.6f}")

print("\n" + "="*60)
print("RECOMMENDATION")
print("="*60)
print("\nFor seismic deconvolution, 'spectral_flatness' is usually best because:")
print("  - It finds wavelets with the broadest, flattest frequency content")
print("  - These wavelets are most stable for inverse filtering")
print("  - Reduces risk of amplifying noise at specific frequencies")
print("\nFor impulsive sources (like hammer hits), 'kurtosis' can work well.")
print("For strong signals, 'energy' ensures you're using the strongest part.")
print("For broadband sources, 'bandwidth' maximizes frequency range.")

print("\n" + "="*60)
print("NEXT STEPS")
print("="*60)
print("\n1. Review the generated diagnostic plots")
print("2. Choose your preferred criterion")
print("3. Update process_xstar.py to use auto_optimize=True:")
print("\n   seis.signature_deconvolution(")
print("       trace_number=6415,")
print("       start_time_ms=30.5,")
print("       end_time_ms=31.7,")
print("       method='wiener',")
print("       epsilon=0.01,")
print("       auto_optimize=True,  # <-- Enable optimization")
print("       optimize_criterion='spectral_flatness')")
print("\n")
