@echo off
REM Double-click (Windows) to build the example analyses and open them.
REM Needs Python (python.org; tick "Add Python to PATH").
cd /d "%~dp0"
for %%C in (pairwise_advanced sglt2_hf nma_inconsistency dta cnma rareevents multivariate rmst) do (
  python run.py configs/example_%%C.json || goto :err
)
echo.
echo Opening reports...
for %%C in (pairwise_advanced sglt2_hf nma_inconsistency dta cnma rareevents multivariate rmst) do (
  start "" "output\example_%%C.html"
)
echo.
echo AACT example (needs a local AACT snapshot; set TRIALFORGE_AACT first):
echo   python run.py configs\example_aact_finerenone.json
echo.
echo Make your own: copy a file in configs\, edit it, then:
echo   python run.py configs\your_file.json
pause
exit /b 0
:err
echo.
echo Could not run. Is Python installed? https://www.python.org/downloads/
pause
exit /b 1
