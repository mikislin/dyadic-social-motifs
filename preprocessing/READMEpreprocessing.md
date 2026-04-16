# Dyadic preprocessing

This package preprocesses dyadic SLEAP tracks with:
- confidence-aware masking from `SLEAPscores`
- 2D jump detection with absolute and quantile thresholds
- short-gap interpolation
- body-geometry repair by marking offending nodes as `NaN` and re-interpolating
- session-specific arena masking
- compact QC outputs using sparse masks
- export and cluster scripts for one-session-per-file preprocessing

## Expected raw input per session
A session struct should contain:
- `SLEAPtracks`: `[T x nodes x 2 x animals]`
- `SLEAPscores`: `[T x nodes]` or `[]`
- `time`: `[T x 1]`
- optional `excludedFrames`: `[T x 1]`

Low-confidence masking is applied as:
```matlab
tracks(scores < 0.2,:,:) = NaN
```
per node and frame.

## Main entry point
```matlab
P = default_preprocessing_params();
out = preprocess_session(sessionRaw, P);
```

## Export dbase to per-session MAT files
```matlab
load('your_dbase.mat', 'dbase')
export_dbase_sessions_to_mat(dbase, '/path/to/preproc_input')
```

## Run preprocessing on exported session files
```matlab
P = default_preprocessing_params();
run_preprocessing_batch('/path/to/preproc_input', '/path/to/preproc_output', P)
```

## Load outputs back into dbase
```matlab
load('your_dbase.mat', 'dbase')
dbase = load_preprocessing_outputs_to_dbase(dbase, '/path/to/preproc_output');
save('your_dbase_preprocessed.mat', 'dbase', '-v7.3')
```

## Cluster deployment
Use `cluster/preprocess_sessions.slurm` as a template.
Workflow:
1. export one MAT file per session
2. submit SLURM array job
3. each task processes one session file
4. merge results back into `dbase`

## Tests
From MATLAB:
```matlab
cd tests
results = runtests('test_preprocess_session.m')
table(results)
```

## Notes
- plotting is off by default for cluster runs
- GPU is not required
- parallelize across sessions, not within a session
