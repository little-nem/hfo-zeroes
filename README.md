Archive structure:

data/
    Contains some marked EEG segments. Those segments come from the archive
    http://meg.univ-amu.fr/data_papers/roehri2017/simulation_rate_3.zip. Their
    simulation process is described in [1]. They are exported to Matlab .mat,
    files using the AnyWave software (http://meg.univ-amu.fr/wiki/AnyWave).

    Note that on Linux, to run AnyWave from its folder, I need to use:

    `env QT_PLUGIN_PATH=plugins/ ./AnyWaveLinux`

demo_matlab/ and demo_octave/
    Contain the “proof of concept” detector described in the report,
    respectively optimised for Matlab and for Octave (the main differences
    are in the .mex and .oct files), the scripts are mostly the same thanks
    to Octave being compatible with the syntax of Matlab.
    From a Matlab or an Octave instance having the Time-Frequency Toolbox
    (http://tftb.nongnu.org/) avalaible on its path, just run:

    ```
    >> detector_demo
    ```

stats_matlab/ and stats_octave/
    Contain a script that produces statistics on a number of marked HFOs.
    Same remark as for demo_{matlab,octave} regarding the diffrences between
    the two folders.
    Here run:

    ```
    >> detector_stats
    ```

# Native Files (mex / oct)

If needs be, one can recompile the native code using either: 

from Matlab:
```
>> mex CXXFLAGS='\$CXXFLAGS -std=c++11" find_components_mex.cpp
>> mex find_zeros_mex.cpp
```

of from Octave:
```
> mkoctfile find_components_oct.cpp
> mkoctfile find_zeros_oct.cpp
```

Refs:
[1] : Roehri, N., Pizzo, F., Bartolomei, F., Wendling, F., & Bénar, C. G. (2017)
    What are the assets and weaknesses of HFO detectors? A benchmark framework
    based on realistic simulations.
    PloS one, 12(4), e0174702.
