🦠 Este repositorio contiene la implementación en MATLAB del modelo SIR con vacunación para el estudio numérico de la propagación de un virus en una población de 7000 individuos. El modelo considera tres compartimentos: susceptibles, infectados y recuperados, e incorpora una tasa de vacunación gamma.

⚙️ Se implementan y comparan tres métodos numéricos: Euler explícito, Runge-Kutta de cuarto orden y la función adaptativa ode45 de MATLAB, todos aplicados sobre el intervalo [0,10] con 500 subintervalos y paso h=0.02.

📁 El script principal Laboratorio_Grupal.m contiene todo el código organizado en secciones comentadas. Las funciones euler_vector y rk4 se incluyen como funciones locales al final del mismo archivo. Los parámetros del modelo son beta=0.0022, mu=0.0137, theta=0.4477 y gamma=0.04, con condiciones iniciales S(0)=6999, I(0)=1 y R(0)=0.

📊 Los resultados muestran que RK4 reproduce prácticamente la solución de referencia de ode45 con errores del orden de 1e-4, mientras que Euler presenta errores del orden de 1e-1. El análisis del efecto de la tasa de vacunación muestra que aumentar gamma de 0.04 a 0.8 reduce el pico de infectados de 5914 a 3161 individuos, una disminución del 47%.

🎓 Este trabajo corresponde al Laboratorio Grupal de la asignatura Métodos Numéricos Aplicados del Máster Universitario en Ingeniería Matemática y Computación de UNIR.
