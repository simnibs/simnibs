set INSTDIR=%1
call %INSTDIR%\simnibs_env\Scripts\activate
python %INSTDIR%\postinstall_simnibs.py -d %INSTDIR% --copy-matlab --setup-links --no-extra-coils