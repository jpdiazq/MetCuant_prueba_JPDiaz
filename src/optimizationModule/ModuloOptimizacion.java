package optimizationModule;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import org.jblas.DoubleMatrix;

public class ModuloOptimizacion {
	
    // -----------------------------------------------------------------
    // Atributos
    // -----------------------------------------------------------------
    
    /**
     * Número de activos.
     */
    
    private int numActivos;
    
    /**
     * Número de datos (Filas).
     */
    
    private int numDatos;
    
    /**
     * Lista Activos.
     */
    
    private Activo[] listaActivos;
    
    // -----------------------------------------------------------------
    // Constructores
    // -----------------------------------------------------------------
        
	/**
	 * Constructor de la clase ModuloOptimizacion. Lee el archivo de CSV que debe estar en la misma carpeta del proyecto y lo guarda en los atributos. Inicializa los demás atributos. 
	 * @throws Exception si existe algún error leyendo el CSV
	 */
	public ModuloOptimizacion(  ) throws Exception
	{		        
		String line = "";  
		String splitBy = ",";
		boolean primeraLinea = true;
		double[][] tempRetornos = new double[1][1];
		String[] tempTickers = new String[1];
				
		try   
		{  
			BufferedReader br = new BufferedReader(new FileReader("retornos_acciones.csv")); 
			
			//contar filas
			int alto = -1;
			int ancho = -1;
			while ((line = br.readLine()) != null)   //returns a Boolean value  
			{  
				alto++; 
			}
						
			br = new BufferedReader(new FileReader("retornos_acciones.csv"));
			
			//Guardar los datos de cada fila
			int contFilas = 0;
			while ((line = br.readLine()) != null)  
			{  
				String[] datos = line.split(splitBy);
				if (primeraLinea == true) {
					ancho = datos.length-1;
					tempTickers = new String[ancho];
					for(int i = 1 ; i<ancho+1 ; i++) {
						tempTickers[i-1] = datos[i];
					}
					tempRetornos = new double[alto][ancho];
					primeraLinea = false;
				}
				else {
					for(int i = 1 ; i < datos.length; i++) {
						tempRetornos[contFilas][i-1] = Double.parseDouble(datos[i]);
					}
					contFilas++;					
				} 
			}  
			
			numActivos = ancho;
			numDatos = alto;
			listaActivos = new Activo[numActivos];
			
			//inicializar activos
			for(int i = 0 ; i<ancho; i++) {
				String tempTicker = tempTickers[i];
				double[] tempRetornosTicker = new double[alto];
				for(int j =0 ; j < alto; j++) {
					tempRetornosTicker[j] = tempRetornos[j][i];
				}
				listaActivos[i] = new Activo(tempTicker, tempRetornosTicker);
			}		
			
			br.close();
		}   
		catch (IOException e)   
		{  
			e.printStackTrace();  
		} 
	}

	
	// -----------------------------------------------------------------
    // Métodos
    // -----------------------------------------------------------------
    
	/**
	 * Retorna el atributo de la lista de activos
	 * @return la lista de activos
	 */
	public Activo[] darListaActivos() {
		return listaActivos;
	}
	
	/**
	 * Calcula la matriz de covarianza para la lista de activos (atributo), dado un inicio del periodo de tiempo por parámetro y una ventana de tiempo establecida.
	 * @param iniInd Inicio del periodo de tiempo para calcular la matriz de covarianza (índice de la serie de tiempo)
	 * @param ventana Ventana de datos para calcular la matriz de covarianza (en número de periodos)
	 * @return la matriz de covarianza calculada para la ventana de tiempo indicada.
	 */
	public double[][] calcularMatrizCov(int iniInd, int ventana) {
		double[][] matrizCov = new double[numActivos][numActivos];
		int finInd = Math.min(iniInd + ventana,numDatos);
		
		try {
			for(int i = 0; i < numActivos; i++) {
				for(int j = 0; j < numActivos; j++) {
					double suma = 0;
					for(int k = iniInd; k < finInd; k++) {
						suma+= (listaActivos[i].darListaRetornos()[k]-listaActivos[i].calcularRetornoMedioVentana(iniInd, ventana))*(listaActivos[j].darListaRetornos()[k]-listaActivos[j].calcularRetornoMedioVentana(iniInd, ventana));
					}
					matrizCov[i][j] = suma/(ventana-1);
				}
			}
		}
		catch (Exception e) {			
			e.printStackTrace();
		}
		
		return matrizCov;
	}
	
	/**
	 * Extrae los retornos para la lista de activos (atributo), dado un inicio del periodo de tiempo por parámetro y una ventana de tiempo establecida.
	 * @param iniInd Inicio del periodo de tiempo para calcular la matriz de retornos (índice de la serie de tiempo)
	 * @param ventana Ventana de datos para calcular la matriz de retornos (en número de periodos)
	 * @return la matriz de retornos extraída para la ventana de tiempo indicada.
	 */
	public double[][] retornosVentana(int iniInd, int ventana) {
		double[][] retornosV = new double[numActivos][1];
		
		try {
			for(int i = 0; i < numActivos; i++) {
				retornosV[i][0] =  listaActivos[i].calcularRetornoMedioVentana(iniInd, ventana);
			}
		}
		catch (Exception e) {			
			e.printStackTrace();
		}
		
		return retornosV;
	}
	
	
	
	/**
	 * Calcula la matriz Q que será utilizada para el método del gradiente conjugado.
	 * @param iniInd Inicio del periodo de tiempo para calcular la matriz Q (índice de la serie de tiempo)
	 * @param ventana Ventana de datos para calcular la matriz Q (en número de periodos)
	 * @return la matriz Q calculada.
	 * @throws Exception si existe algún error al calcular el retorno medio según una ventana de tiempo (que el índice de inicio sea mayor al número de datos.
	 */
	public DoubleMatrix mat_Q(int iniInd, int ventana) throws Exception {
		
		double[][] matrizCov = calcularMatrizCov(iniInd, ventana);
		
		int size_Q = matrizCov.length + 2;
		
		double[][] matriz_Q = new double[size_Q][size_Q];
		
		for(int i = 0; i<size_Q; i++) {
			for(int j = 0; j<size_Q; j++) {
				if(i<matrizCov.length && j<matrizCov.length) {
					matriz_Q[i][j] = matrizCov[i][j];
				}
				else if(i<matrizCov.length && j==matrizCov.length) {
					matriz_Q[i][j] = listaActivos[i].calcularRetornoMedioVentana(iniInd, ventana);
				}
				else if(i<matrizCov.length && j==(size_Q-1)) {
					matriz_Q[i][j] = 1;
				}
				else if(i==matrizCov.length && j<matrizCov.length) {
					matriz_Q[i][j] = listaActivos[j].calcularRetornoMedioVentana(iniInd, ventana);
				}
				else if(i==(size_Q-1) && j<matrizCov.length) {
					matriz_Q[i][j] = 1;
				}
				else {
					matriz_Q[i][j] = 0;
				}				
			}
		}
		
		return new DoubleMatrix(matriz_Q);		
	}
	
	/**
	 * Calcula el vector x_0 para el método de gradiente conjugado
	 * @param x_inicial Escalar que determina el valor para inicializar x en el método de gradiente conjugado
	 * @return Vector completo x_0 para inicializar x en el método de gradiente conjugado
	 */
	public DoubleMatrix vect_x0(double x_inicial){
			
			double[][] vectx0 = new double[listaActivos.length+2][1];
			
			for(int i = 0; i< listaActivos.length+2; i++) {
				vectx0[i][0] = x_inicial;
			}
			
			return new DoubleMatrix(vectx0);		
	}
	
	/**
	 * Calcula el vector b para el método de gradiente conjugado
	 * @param r_target Escalar que determina el retorno objetivo.
	 * @return Vector completo b para el método de gradiente conjugado
	 */
	public DoubleMatrix vect_b(double r_target) {
		
		double[][] vectb = new double[listaActivos.length+2][1];
		
		for(int i = 0; i< listaActivos.length; i++) {
			vectb[i][0] = 0;			
		}
		
		vectb[listaActivos.length][0] = r_target;
		vectb[listaActivos.length+1][0] = 1;
		
		return new DoubleMatrix(vectb);		
}
	
	/**
	 * Método para multiplicar una matriz por un escalar.
	 * @param S matriz a ser multiplicada
	 * @param scalar escalar
	 * @return Matriz multiplicada por el escalar.
	 */
	public DoubleMatrix mat_mulScal(DoubleMatrix S, double scalar) {
		
		double[][] ret_mat = new double[S.rows][S.columns];
		
		for (int i=0; i<S.rows; i++) {
			
			for(int j=0; j<S.columns; j++) {
				
				ret_mat[i][j] = S.get(i,j)*scalar;				
			}
		}
		
		return new DoubleMatrix(ret_mat);
		
	}
	
	/**
	 * Función que realiza el algoritmo del gradiente conjugado
	 * @param Q Matriz Q del sistema lineal de ecuaciones
	 * @param b Vector b del sistema lineal de ecuaciones
	 * @param x_0 x inicial del sistema lineal de ecuaciones
	 * @param tol Tolerancia para detener el algoritmo
	 * @param maxIter Máximo número de iteraciones, también detiene el algoritmo
	 * @return El X resultado del algoritmo del gradiente conjugado.
	 */
	public double[][] gradiente_conjugado(DoubleMatrix Q, DoubleMatrix b, DoubleMatrix x_0, double tol, int maxIter) {
		
		double[][] x_mat = new double[Q.rows][1];
		
		//Inicialización
		DoubleMatrix x = x_0;
		DoubleMatrix s = b.sub(Q.mmul(x_0));
		DoubleMatrix p = b.sub(Q.mmul(x_0));
		DoubleMatrix s_ant = s.transpose().mmul(s);
		double residual_old = s_ant.get(0,0);		
		
		//Iteraciones
		int iter = 0;
		
		while(iter<maxIter) {
			
			//Alpha
			DoubleMatrix Q_p = Q.mmul(p);
			double alpha = residual_old/p.transpose().mmul(Q_p).get(0,0);
						
			//X
			x = x.add(mat_mulScal(p,alpha));
						
			//S
			s = s.sub(mat_mulScal(Q_p,alpha));
			
			//residual_num
			double residual_new = s.transpose().mmul(s).get(0,0);
			
			for (int i=0; i<x.rows;i++) {
				x_mat[i][0] = x.get(i,0);
			}
			
			//Revisar condición de parada
			if (Math.sqrt(residual_new)<=tol) {
				break;
			}
						
			//P
			p = s.add(mat_mulScal(p,(residual_new/residual_old))); //Beta
			
			//Actualizar residual viejo
			residual_old = residual_new;
			
			iter++;
		}
				
		return x_mat;
				
	}
	
	//Agregar retorno arbitrario y tolerancia como parámetro
	/**
	 * Realiza la optimización de un portafolio para un periodo de tiempo, establecido por un índice inicial y una ventana de tiempo. También está sujeto por un retorno objetivo.
	 * @param iniInd Indice que establece la posición inicial del periodo de tiempo para la optimización
	 * @param ventana Ventana de teimpo para la optimización
	 * @param retObjetivo Retorno objetivo para la optimización (minimización de varianza)
	 * @return La solución de la optimización x_sol, dado un periodo de tiempo
	 * @throws Exception En caso de que el índice de inicio sea mayor al número de datos
	 */
	public double[][] optimizacionPortafolio(int iniInd, int ventana, double retObjetivo) throws Exception {
		
		
		double tolerancia = 0.00000001;
		double x_inicial = 0.0001;
		int maxIter = 10000;
		
		DoubleMatrix Q = mat_Q(iniInd, ventana);
		DoubleMatrix x0 = vect_x0(x_inicial);
		DoubleMatrix b = vect_b(retObjetivo);
		
		double[][] x_solucion = gradiente_conjugado(Q, b, x0, tolerancia, maxIter);
		
		return x_solucion;
	}
	
	//Obtener sólo pesos del portafolio
	/**
	 * Extrae sólo los pesos 'w' del vector de solución x (el cual también contiene los multiplicadores de Lagrange)
	 * @param x_completo El vector completo, resultado de la optimización de portafolio.
	 * @return Los pesos 'w'
	 */
	public DoubleMatrix soloPesos(DoubleMatrix x_completo) {
		
		double[][] soloPesosRet = new double[x_completo.rows-2][1];
		
		for (int i=0; i<x_completo.rows-2; i++) {	
			soloPesosRet[i][0] = x_completo.get(i,0);
		}
		
		return new DoubleMatrix(soloPesosRet);
		
	}
	
	
	//Calcular Covarianza
	/**
	 * Calcula la covarianza de un portafolio para un vector de pesos 'w' y una matriz de covarianza 'Sigma'
	 * @param covMat Matriz de covarianza del portafolio
	 * @param pesos Pesos para calcular la covarianza del portafolio
	 * @return Covarianza del portafolio.
	 */
	public double calcCovar(DoubleMatrix covMat, DoubleMatrix pesos) {
		double respCalcCovar = 0;
		
		respCalcCovar = pesos.transpose().mmul(covMat.mmul(pesos)).get(0,0);
		
		return respCalcCovar;
	}
	
	
	/**
	 * Método para imprimir matrices del tipo DoubleMatrix
	 * @param mat Matrix para ser impresa en la consola
	 */
	public void print_mat(DoubleMatrix mat) {
		//For traversing through rows
		for(int i =0;i<mat.rows;i++) {
			//For traversing through columns
			for(int j=0;j<mat.columns;j++) {
				
				System.out.print("\t"+mat.get(i,j)+"\t");
				
			}
			System.out.println();
		}
		
	}
	

	
	//Backtesting
	
	/**
	 * Método que realiza el backtesting de la estrategia de optimización, dada una ventana de tiempo 'InSample' y una ventana de tiempo 'OutSample'
	 * @param nInSample Número de periodos para utilizar como conjunto 'In Sample'
	 * @param nOutSample Número de periodos para utilizar como conjunto 'Out Sample'
	 * @param nStep Pasos para avanzar en el backtesting (periodos que se mueven ambas ventanas)
	 * @param dailyRfRate Tasa libre de riesgo con la misma periodicidad de los datos, para calcular Sharpe Ratio
	 * @param resultsOutSample Indica si se quieren obtener los resultados promedio del backtesting para la ventana 'InSample' (false) o 'OutSample' (true)
	 * @return Una matriz con el resultado del backtesting, que obtiene los promedios de la ventana elegida para el resultado, de las siguientes medidas: 'Target Return' 'Average Return' 'Average Variance' 'Average Sharpe Ratio'
	 * @throws Exception En caso de que el índice de inicio sea mayor al número de datos
	 */
	public double[][] backtesting(int nInSample, int nOutSample, int nStep, double dailyRfRate, boolean resultsOutSample) throws Exception{
		double[] retObjetivo = {0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.005,0.01,0.02,0.05,0.1};
		int posibIter = (int) Math.ceil((numDatos - (nInSample+nOutSample))/nOutSample);

		
		double[][] resultsReturnSample = new double[retObjetivo.length][posibIter];
		double[][] resultsVolSample = new double[retObjetivo.length][posibIter];
		double[][] resultsSharpeSample = new double[retObjetivo.length][posibIter];
		
		int iter = 0;
		
		int posIniInSample = 0;
		int posIniOutSample = posIniInSample+nInSample;
		
		while (iter < posibIter) {
			
			for(int i=0; i<retObjetivo.length; i++) {
				//Optimización en Periodo IN Sample
				double[][] x_sol = optimizacionPortafolio(posIniInSample, nInSample,retObjetivo[i]);
				DoubleMatrix x_mat = new DoubleMatrix(x_sol);
				
				DoubleMatrix retVentIn =  new DoubleMatrix(retornosVentana(posIniInSample, nInSample));
				DoubleMatrix covMatIn = new DoubleMatrix(calcularMatrizCov(posIniInSample, nInSample));
				double tempReturnIn = retVentIn.transpose().mmul(soloPesos(x_mat)).get(0,0);
				double tempVolIn = Math.sqrt(calcCovar(covMatIn, soloPesos(x_mat)));
				double tempSharpeIn = (tempReturnIn-dailyRfRate)/tempVolIn;
				
				//Cálculos Periodo Out Sample
				
				DoubleMatrix retVent =  new DoubleMatrix(retornosVentana(posIniOutSample, nOutSample));
				DoubleMatrix covMat = new DoubleMatrix(calcularMatrizCov(posIniOutSample, nOutSample));
				double tempReturnOut = retVent.transpose().mmul(soloPesos(x_mat)).get(0,0);
				double tempVolOut =  Math.sqrt(calcCovar(covMat, soloPesos(x_mat)));
				double tempSharpeOut = (tempReturnOut-dailyRfRate)/tempVolOut;
				
				if(resultsOutSample) {
					resultsReturnSample[i][iter] = tempReturnOut;
					resultsVolSample[i][iter] = tempVolIn;
					resultsSharpeSample[i][iter] = tempSharpeOut;
				}
				else {
					resultsReturnSample[i][iter] = tempReturnIn;
					resultsVolSample[i][iter] = tempVolOut;
					resultsSharpeSample[i][iter] = tempSharpeIn;
				}
				
			}
			
			posIniInSample += nStep;
			posIniOutSample = posIniInSample+nInSample;
			iter++;
		}
		
		double[][] avgResults = new double[retObjetivo.length][4];
		
		for(int i = 0; i<retObjetivo.length; i++) {
			avgResults[i][0] = retObjetivo[i];
			for(int j = 1; j<4; j++) {
				double[][] arrayToUse = null;
				switch(j) {
				  case 1:
					arrayToUse = resultsReturnSample;
				    break;
				  case 2:
				    arrayToUse = resultsVolSample;
					break;
				  case 3:
					arrayToUse = resultsSharpeSample;
				    break;
				}
				
				double[] subSetSelected = new double[posibIter];
				for(int k=0; k<posibIter; k++) {
					subSetSelected[k] = arrayToUse[i][k];
				}
				
				double avgMeasure = average(subSetSelected);
				avgResults[i][j] = avgMeasure;				
			}
		}
		
		return avgResults;
		
	}
	
	/**
	 * Calcula el promedio para un arreglo de datos
	 * @param arreglo Vector para calcular el promedio
	 * @return Promedio del vector pasado por parámetro
	 */
	public double average(double[] arreglo) {
		double suma = 0;
		for(int i = 0; i<arreglo.length; i++) {
			suma+=arreglo[i];
		}
		return suma/arreglo.length;
	}
	
	// -----------------------------------------------------------------
    // Main
    // -----------------------------------------------------------------
    
	/**
	 * Método main
	 * @param args
	 */
	public static void main(String[] args) {
		
		try {
			ModuloOptimizacion opti = new ModuloOptimizacion();
			
			System.out.println("Corriendo backtesting Out Sample...");
			double[][] avgResOutSample = opti.backtesting(100, 12, 12, 0.00005,true);
			
			System.out.println("Corriendo backtesting In Sample...");
			double[][] avgResInSample = opti.backtesting(100, 12, 12, 0.00005,false);
			
			DoubleMatrix avgResDMOut = new DoubleMatrix(avgResOutSample);
			DoubleMatrix avgResDMIn = new DoubleMatrix(avgResInSample);
									
			System.out.println("Avg Results Out Sample (Tgt Return, Avg Return, Avg Variance and Avg Sharpe Ratio):");
			opti.print_mat(avgResDMOut);	
			
			
			System.out.println("Avg Results In Sample (Tgt Return, Avg Return, Avg Variance and Avg Sharpe Ratio):");
			opti.print_mat(avgResDMIn);					
			
		} 
		catch (Exception e) {
			
			e.printStackTrace();
		}
	}

}
