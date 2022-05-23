/**
 * 
 */
package optimizationModule;


/**
 * @author jpdia
 *
 */
public class Activo {

    
    // -----------------------------------------------------------------
    // Atributos
    // -----------------------------------------------------------------
    
    /**
     * Ticker.
     */
    
    private String ticker;
    
    /**
     * Lista de Retornos.
     */
    
    private double[] listaRetornos;
    
    
    // -----------------------------------------------------------------
    // Constructores
    // -----------------------------------------------------------------
    

    
    /**
     * @param paramTicker
     * @param paramRetornos
     * @throws Exception
     */
    public Activo(String paramTicker, double[] paramRetornos ) throws Exception
    {
    	
    	ticker = paramTicker;
    	listaRetornos = paramRetornos;
    	
    }
    
    
    //-----------------------------------------------------------------
    // M�todos
    //-----------------------------------------------------------------
    
    
    
    /**
     * Retorna el ticker
     * @return String ticker del activo.
     */
    
    public String darTicker()
    {
    	return ticker;
    }
    
    /**
     * Retorna la lista de retornos
     * @return double[] lista de retornos hist�ricos del activo.
     */
    
    public double[] darListaRetornos()
    {
    	return listaRetornos;
    }


    
    /**
     * Calcula el retorno medio del activo, dado un periodo de tiempo establecido por un �ndice de inicio y una ventana de datos.
     * @param iniInd �ndice que establece el comienzo del periodo de tiempo para calcular el retorno medio.
     * @param ventana ventana de tiempo para calcular el retorno medio.
     * @return El retorno medio del activo para la ventana indicada.
     * @throws Exception si el �ndice es superior a la cantidad de datos.
     */
    public double calcularRetornoMedioVentana(int iniInd, int ventana) throws Exception
    {
    	double suma = 0;
    	double retornoMedio = 0;
    	if(iniInd >= listaRetornos.length)
    	{
    		throw new Exception("Indice es superior a la cantidad de datos");
    	}
    	else
    	{
    		int finInd = Math.min(iniInd + ventana,listaRetornos.length);
    		for(int i = iniInd; i < finInd; i++ )
            {
        		suma+=listaRetornos[i];
            }
        	retornoMedio = suma/(finInd-iniInd);
    	}
    	
    	return retornoMedio;
    	
    }
    
    
    
    
    
    
    

}
