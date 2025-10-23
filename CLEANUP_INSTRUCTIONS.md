# Cleaning Up Git History - Large Files Removal

## Problem
The `.git` directory is currently ~52 MB, making clones slow. This is caused by large binary files that were committed to the repository.

## Files Updated
- `.gitignore` - Updated to exclude large binary files going forward
- `deltaseis/data/README.md` - Created to document data file handling
- This document - Instructions for cleanup

## Solution

The large binary files have been removed from the working directory and `.gitignore` has been updated to prevent future commits. However, these files still exist in Git history.

### For Repository Maintainers: Cleaning Git History

⚠️ **WARNING**: The following steps rewrite Git history. All contributors must re-clone the repository after this is done.

#### Option 1: Using git-filter-repo (Recommended)

```bash
# Install git-filter-repo
pip install git-filter-repo

# Create a backup first!
git clone --mirror https://github.com/Deltares-research/DeltaSEIS.git DeltaSEIS-backup

# Remove specific file extensions from history
cd DeltaSEIS
git filter-repo --invert-paths --path-glob '*.sgy' --force
git filter-repo --invert-paths --path-glob '*.png' --force
git filter-repo --invert-paths --path-glob '*.avi' --force
git filter-repo --invert-paths --path-glob '*.pkl' --force
git filter-repo --invert-paths --path-glob '*.asc' --force

# Or remove entire directories
git filter-repo --invert-paths --path 'deltaseis/data/' --force
git filter-repo --invert-paths --path 'examples/*.avi' --force

# Force push to update remote (requires admin privileges)
git push origin --force --all
git push origin --force --tags
```

#### Option 2: Using BFG Repo-Cleaner

```bash
# Download BFG
# https://rtyley.github.io/bfg-repo-cleaner/

# Clone a fresh bare repository
git clone --mirror https://github.com/Deltares-research/DeltaSEIS.git

# Remove files larger than 1MB
java -jar bfg.jar --strip-blobs-bigger-than 1M DeltaSEIS.git

# Clean up
cd DeltaSEIS.git
git reflog expire --expire=now --all
git gc --prune=now --aggressive

# Force push (requires admin privileges)
git push --force
```

### After History Cleanup

1. **Verify the size reduction**:
   ```bash
   git count-objects -vH
   ```

2. **Notify all contributors** to:
   - Delete their local clones
   - Re-clone the repository
   - Do NOT merge old branches - rebase or recreate them

3. **Protect main branch** to prevent accidental large file commits

### Expected Results

After cleanup:
- `.git` directory should be < 5 MB (down from 52 MB)
- Clone times should be significantly faster
- All code and small text files remain intact
- Git history is rewritten (new commit SHAs)

## Alternative: Keep History and Use Git LFS

If you want to keep the data files but reduce clone times:

```bash
# Install Git LFS
git lfs install

# Track large file types
git lfs track "*.sgy"
git lfs track "*.png"
git lfs track "*.avi"
git lfs track "*.pkl"

# This will only affect new commits
# You'd still need to clean history as above
```

## Questions?

Contact the repository maintainers if you have questions about this cleanup process.
