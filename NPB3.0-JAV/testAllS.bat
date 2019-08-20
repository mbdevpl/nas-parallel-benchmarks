set CLASSPATH=d:\\jdk1.1.8\\src
set NPBPATH=d:\\cygwin\\home\\frumkin\\java\\ngb


start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.BT -serial CLASS=S


start /b /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.BT -np2 CLASS=S



start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.SP -serial CLASS=S


start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.SP -np2 CLASS=S



start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.LU -serial CLASS=S


start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.LU -np2 CLASS=S



start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.FT -serial CLASS=S


start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.FT -np2 CLASS=S



start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.MG -serial CLASS=S


start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.MG -2 CLASS=S



start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.CG -serial CLASS=S


start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.CG -np2 CLASS=S



start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.IS -serial CLASS=S


start /B /wait java -classpath %NPBPATH% -mx200MB NPB3_0_JAV.IS -np2 CLASS=S


