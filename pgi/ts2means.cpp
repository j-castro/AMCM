/*
Parallel C implementation of the ts2means algorithm, originally
implemented in MATLAB by Professor Arturo Camacho.

JORGE CASTRO - A71602 - UCR

 OBSERVACIONES:
 1- EL TAMANO MAXIMO DE LA VENTANA ES EL TAMANO DE LA SENAL
 2- la parte diferente con respecto al algoritmo en matalab resulta en
 el caso especial en que para la interpolacion se necesita una ventana extra
 al final, resulta en una variacion promedio de un 1% de los datos.
*/


#include "ts2means.h"


//Algoritmo de verificacion de parametros del ts2means, retorna un 1 si los parametros estan bien, sino un 0
int verificacion(double dt, double wincrement, double woverlap, char* BIAS,intervalo wdurationlim) {

    int resultado = 1;

    if (dt <= 0){
        resultado = 0;
        fprintf(stderr, "Time interval must be positive. \n");
    }
    else if (wincrement <= 1){
        resultado = 0;
        fprintf(stderr, "Window increment factor must be greater than 1. \n");
    }
    else if (woverlap>1 || woverlap<0){
        resultado = 0;
        fprintf(stderr, "Window overlap must be between 0 and 1. \n");
    }
    else if (strcmp(BIAS,"max") != 0 && strcmp(BIAS,"min") != 0) {
        resultado = 0;
        fprintf(stderr, "Invalid bias argument. Must be: 'max' or 'min'. \n");
    }
    else if (wdurationlim.inicio > wdurationlim.fin){
        resultado = 0;
        fprintf(stderr, "Invalid window duration limit. \n");
    }

    return resultado;

}


//TS2MEANS SECUENCIAL
//RECIBE:   1. Vector con la senal a clasificar 'x', 2. Intervalo de muestreo 'dt'  3. Intervalo de ventanas validas
//          4. Incremento del tamano de las ventanas 5. Traslape porcentual de las ventanas 6. Para el post-procesamiento
//REALIZA:  el algoritmo de clasificacion en 2 clases de una senal (umbralizacion dinamica)
//RETORNA:  un vector con la senal clasificada
vec ts2means(vec x,double dt,intervalo wdurationlim,double wincrement,double woverlap,char* BIAS){

	vec signal;

    int i,j,m;          //indices para uso de loops
    int indice;         //indice para la distribucion de ventanas en la senal
    int w4signal;       //number of windows of the current size needed for the hole signal
    int i1,i2;          //indices para el ajuste de las ventanas a la senal
    int hopSz;          //tamanho del brinco entre ventanas, expresado en cantidad de muestras
    int scFlag;         //Bndera del caso especial

    vec t;           //vector del tiempo
    vec wsize;       //vector de tamanho de las ventanas en cantidad de muestras
    vec wduration;   //vector de duracion de las ventanas en segundos
    vec ti;          //time points to center de hann window
    vec xnp;         //vector with NaN's, x and more NaN's
    vec hannW;       //the hann window
    vec tempHW;      //temporary hann window
    vec tempX;       //temporary X fragment
    vec temp1, temp2, temp3;	    //vectores temporales miscelaneos
    vec ci0;         //centroids class 0  of and specific window over all X
    vec ci1;         //centroids class 1  of and specific window over all X

    matrix D;           //matriz de distancias
    matrix C0;          //matriz de centroides clase 0
    matrix C1;          //matriz de centroides clase 1

    //Verificacion de que los parametros esten correctos
    if (verificacion(dt,wincrement,woverlap,BIAS,wdurationlim) == 1) {

        //printf("parametros bien ajustados.  dt = %f, winc = %f, wov = %f, BIAS : %s \n",dt,wincrement,woverlap,BIAS);

        //#############################################################################################################
        //PARTE 1 - OBTENCION DEL VECTOR DE TAMANHO DE VENTANAS Y DURACION DE LAS MISMAS
        wsize = getWSvec(wdurationlim, wincrement, dt, x.x);
        wduration = makev(wsize.x);

        for(i = 0;i < wsize.x ;i++){
            wduration.v[i] = wsize.v[i] * dt; //vector de duracion de las ventanas
        }

        //printf("El vector de tamanho de ventanas en cantidad de muestras es: \n");
        //printv(wsize);

        //printf("El vector de tamanho de ventanas en cantidad de segundos es: \n");
        //printv(wduration);

        //FIN PARTE 1
        //#############################################################################################################

        //#############################################################################################################
        //PARTE 2, OBTENCION DE LA MATRIZ DE CENTROIDES Y DISTANCIAS

        //Creacion de las matrices de centroides y distancias entre los centroides
        t = makeFillv(0,dt,x.x);
        D = zerom(wsize.x, t.x);
        C0 = zerom(wsize.x, t.x);
        C1 = zerom( wsize.x, t.x);

        //Calculo de centroides para cada ventana********************************************************************************
        for(i = 0;i < wsize.x; i++) {

            //printf("\n------------CICLO EXTERNO %d\n\n",i);
            scFlag = 0; //bandera para el caso especial

            hopSz = (int) max( 1 , round( wsize.v[i] * (double)(1-woverlap)) );   //tamanho en muestras, del brinco entre ventanas

            w4signal = (int)( (x.x + wsize.v[i]/2 )/hopSz);    //cantidad de ventanas que se ocuparan para cubrir toda la muestra (1 extra for interpolation)
	    //printf("Windows for signal: %d\n",w4signal);
            if ((int)(x.x + wsize.v[i]/2 )% hopSz > 0){w4signal++;}

            ti = makeFillv(0, (double)hopSz*dt, w4signal);   //times in wich hann window is evaluated (middpoint not the start)
                                                            //for later interpolation
            xnp = makev(wsize.v[i]/2 + x.x + wsize.v[i]);   //one extra window at the end to cover all windows positions

            for(m = 0;m < xnp.x; m++){
                if ((m < (wsize.v[i]/2)) || (m >= (wsize.v[i]/2 + x.x))){xnp.v[m] = NAN;} //NaN's X NaN's --> xnp
                else { xnp.v[m] = x.v[ m-(int)(wsize.v[i]/2) ];}
            }

            //ventana hann
            hannW = makeHannv(wsize.v[i]);

            //Verificacion caso especial en el cual ocupo un centroide extra al final de la muestra
            if (ti.v[ti.x - 1] < t.v[t.x - 1]){

                ti = addElementv(ti, t.v[t.x - 1] );  //agrego un ultimo valor para interpolacion
                scFlag = 1;

            }

            //centroides para el tamaño de ventana actual
            ci0 = makev(w4signal + scFlag);
            ci1 = makev(w4signal + scFlag);

            indice = 0; //indice para desplazarme en xnp

            //Desplazamiento de ventanas a través del tiempo según el woverlap!
            for(j = 0;j < w4signal; j++) {

                //printf("------------CICLO INTERNO %d\n",j);

                //filtrado de NaN's ---------------------------------------------------------------------------
                indice = j*hopSz;   //hop multiplied by # of iteration
                i2 = j*hopSz + wsize.v[i];  //the end position of the actual window

                //First NaN zone
                while (isnan(xnp.v[indice])) {indice++;}
                i1 = indice; //first index

                //normal numbers
                while(isnan(xnp.v[indice]) != 1 && indice < i2){ indice++;}

                i2 = indice-1;

		
                tempX = getSubv(xnp,i1,i2);
                tempHW = getSubv(hannW,i1-(j*hopSz),i2-(j*hopSz));
                //--------------------------------------------------------------------------------------------------

                //w2means
                temp1 = w2means(tempX,tempHW);

				ci0.v[j] = temp1.v[0];
				ci1.v[j] = temp1.v[1];

				freev(temp1);
                freev(tempX);
                freev(tempHW);
            }

            //Special Case
            if (scFlag == 1){   //Compute centroids centered at the end of the signal (for later interpolation)
		
                //Calculate the last centroids centered in t.v[t.x -1]
                i1 = x.x - wsize.v[i]/2;
                i2 = x.x - 1;

                tempX = getSubv(x,i1,i2);
                tempHW = getSubv(hannW, 0, (wsize.v[i]/2) -1);

                temp1 = w2means(tempX,tempHW); printf("Centroides obtenidos CASO ESPECIAL: C0: %f  C1: %f\n",temp1.v[0],temp1.v[1]);

				ci0.v[ci0.x -1] = temp1.v[0];
				ci1.v[ci1.x -1] = temp1.v[1];

                //-------liberacion de memoria---------------
                freev(temp1);   freev(tempX);   freev(tempHW);
            }

            //INTERPOLACION LINEAL DE VECTORES********
            temp1 = interp1(ti, ci0, t);
            copyRow(&C0,temp1,i);

            temp2 = interp1(ti, ci1, t);
            copyRow(&C1,temp2,i);

            temp3 = distancev(temp1,temp2);
            copyRow(&D,temp3,i);
            //***************************************

            //liberación de memoria
            freev(temp1);
            freev(temp2);
            freev(temp3);
            freev(ti);
            freev(xnp);
            freev(hannW);
            freev(ci0);
            freev(ci1);

        }
        //***********************************************************************************************************************

        //FIN PARTE 2
        //#############################################################################################################

        //#############################################################################################################
        //PARTE 3 DETERMINACION DEL UMBRAL DE CLASIFICACION Y POST-PROCESAMIENTO

        //Eleccion de los mejores centroides para cada muestra en el tiempo, calculo de su punto medio y umbralizacion
        signal = classifySig(&x,&C0,&C1,&D);
	
        //Post-Procesamiento =)
        postProcesing(&x,&signal,BIAS);

        //Liberar toda la memoria utilizada
        freem(D);
        freem(C0);
        freem(C1);
        freev(t);
        freev(wsize);
        freev(wduration);

    }

	return signal;
}

//TS2MEANS PARALELIZADO
//RECIBE:   1. Vector con la senal a clasificar 'x', 2. Intervalo de muestreo 'dt'  3. Intervalo de ventanas validas
//          4. Incremento del tamano de las ventanas 5. Traslape porcentual de las ventanas 6. Para el post-procesamiento
//REALIZA:  el algoritmo de clasificacion en 2 clases de una senal (umbralizacion dinamica) paralelizado con openmp
//RETORNA:  un vector con la senal clasificada
vec ts2meansP(vec x,double dt,intervalo wdurationlim,double wincrement,double woverlap,char* BIAS){

	vec signal;	    //senial de salida

    int nthreads, tid;      //numero de hilos y su respectivo id
    int i,j,m;              //indices para uso de loops
    int indice;             //indice para la distribucion de ventanas en la senal
    int w4signal;           //number of windows of the current size needed for the hole signal
    int i1,i2;              //indices para el ajuste de las ventanas a la senal
    int hopSz;              //tamanho del brinco entre ventanas, expresado en cantidad de muestras
    int scFlag;             //FLAG for special case
    int iniOffset;	    //Initial offset for processing the first window

    vec t;               //vector del tiempo
    vec wsize;           //vector de tamanho de las ventanas en cantidad de muestras
    vec wduration;       //vector de duracion de las ventanas en segundos
    vec ti;              //time points to center de hann window
    vec xnp;             //vector with NaN's, x and more NaN's
    vec hannW;           //the hann window
    vec tempHW;          //temporary hann window
    vec tempX;           //temporary X fragment
    vec temp1, temp2, temp3;	//vectores temporales miscelaneos
    vec ci0;             //centroides clase 0  de la iesima ventana a lo largo de X
    vec ci1;             //centroids clase 1 of and specific window over all X

    matrix D;               //matriz de distancias
    matrix C0;              //matriz de centroides clase 0
    matrix C1;              //matriz de centroides clase 1


    //Verificacion de que los parametros esten correctos
    if (verificacion(dt,wincrement,woverlap,BIAS,wdurationlim) == 1) {

        //printf("parametros bien ajustados.  dt = %f, winc = %f, wov = %f, BIAS : %s \n",dt,wincrement,woverlap,BIAS);

        //#############################################################################################################
        //PARTE 1 - OBTENCION DEL VECTOR DE TAMANHO DE VENTANAS Y DURACION DE LAS MISMAS

        wsize = getWSvec(wdurationlim, wincrement, dt, x.x);
        wduration = makev(wsize.x);

        for(i = 0;i < wsize.x ;i++){
            wduration.v[i] = wsize.v[i] * dt; //vector de duracion de las ventanas
        }

        //printf("El vector de tamanho de ventanas en cantidad de muestras es: \n");
        //printv(wsize);

        //printf("El vector de tamanho de ventanas en cantidad de segundos es: \n");
        //printv(wduration);

        //FIN PARTE 1
        //#############################################################################################################

        //#############################################################################################################
        //PARTE 2, OBTENCION DE LA MATRIZ DE CENTROIDES Y DISTANCIAS

        //Creacion de las matrices de centroides y distancias entre los centroides
        t = makeFillv(0,dt,x.x);
        D = zerom(wsize.x, t.x);
        C0 = zerom(wsize.x, t.x);
        C1 = zerom( wsize.x, t.x);

	xnp = makev(wsize.v[wsize.x -1]/2 + x.x + wsize.v[wsize.x -1]);   //one extra window at the end to cover all windows positions
	
        for(m = 0;m < xnp.x; m++){
            if ((m < (wsize.v[wsize.x -1]/2)) || (m >= (wsize.v[wsize.x -1]/2 + x.x))){xnp.v[m] = NAN;} //NaN's X NaN's --> xnp
            else { xnp.v[m] = x.v[ m-(int)(wsize.v[wsize.x -1]/2) ];}
        }
        
        double tInicio, tFinal;
        tInicio = omp_get_wtime();
        //#pragma omp parallel private(nthreads, tid, scFlag, hopSz, w4signal, iniOffset, ti, hannW, ci0, ci1, i, j, m, indice, i1, i2, tempX, tempHW, temp1, temp2, temp3) //Paralelizacion

        {

            // Get thread number
            tid = omp_get_thread_num();
            // Get number of threads
            nthreads = omp_get_num_threads();

            /* Only master thread does this*/
            if (tid == 0) {
              printf("Numero de hilos = %d\n", nthreads);
	      printf("tamanho del vector de ventanas wsize = %d\n", wsize.x);
            }

            i = tid;
            //Compute centroids for each window size********************************************************************************
            while(i < wsize.x) {

                //printf("\n------------CICLO EXTERNO %d\n\n",i);
                scFlag = 0; //FLAG for the special case

                hopSz = (int) max( 1 , round( wsize.v[i] * (double)(1-woverlap)) );   //tamanho en muestras, del brinco entre ventanas

                w4signal = (int)( (x.x + wsize.v[i]/2 )/hopSz);    //cantidad de ventanas que se ocuparan para cubrir toda la muestra (1 extra for interpolation)
                if ((int)(x.x + wsize.v[i]/2 )% hopSz > 0){w4signal++;}

                ti = makeFillv(0, (double)hopSz*dt, w4signal);   //times in wich hann window is evaluated (middpoint not the start)
                                                                //for later interpolation

                //ventana hann
                hannW = makeHannv(wsize.v[i]);

                //Verificacion caso especial en el cual ocupo un centroide extra al final de la muestra
                if (ti.v[ti.x - 1] < t.v[t.x - 1]){

                    //printf("Necesito una ventana extra para luego interpolar\n");
                    //system("pause");

                    ti = addElementv(ti, t.v[t.x - 1] );  //agrego un ultimo valor para interpolacion
                    scFlag = 1;

                }

                //centroides para el tamaño de ventana actual
                ci0 = makev(w4signal + scFlag);
                ci1 = makev(w4signal + scFlag);

                indice = 0; //indice para desplazarme en xnp

		//TODO Desplazamiento de las ventanas segun su tamaño 
		iniOffset = wsize.v[wsize.x -1]/2 - wsize.v[i]/2;
		
		//printf("Initial offset = %d of thread %d\n",iniOffset,i);		
		//#pragma omp barrier

                //Desplazamiento de ventanas a través del tiempo según el woverlap!
                for(j = 0;j < w4signal; j++) {

                   // printf("------------CICLO INTERNO %d\n",j);

                    //filtrado de NaN's ---------------------------------------------------------------------------
                    indice = j*hopSz + iniOffset;   //hop multiplied by # of iteration
                    i2 = j*hopSz + iniOffset + wsize.v[i];  //the end position of the actual window

                    //First NaN zone
                    while (isnan(xnp.v[indice])) {indice++;}
                    i1 = indice; //first index

                    //normal numbers
                    while(isnan(xnp.v[indice]) != 1 && indice < i2){ indice++;}

                    i2 = indice-1;
                    tempX = getSubv(xnp,i1,i2);
                    tempHW = getSubv(hannW,i1-(j*hopSz+iniOffset),i2-(j*hopSz+iniOffset));
                    //--------------------------------------------------------------------------------------------------

                    //w2means
                    temp1 = w2means(tempX,tempHW);
                    //printf("Centroides obtenidos: C0: %f  C1: %f\n",temp1.v[0],temp1.v[1]);

                    ci0.v[j] = temp1.v[0];
                    ci1.v[j] = temp1.v[1];

                    freev(temp1);
                    freev(tempX);
                    freev(tempHW);
                }

                //Special Case
                if (scFlag == 1){

                    //Calculate the last centroids centered in t.v[t.x -1]
                    i1 = x.x - wsize.v[i]/2;
                    i2 = x.x - 1;

                    tempX = getSubv(x,i1,i2);
                    tempHW = getSubv(hannW, 0, (wsize.v[i]/2) -1);

                    //w2means
                    temp1 = w2means(tempX,tempHW); printf("Centroides obtenidos CASO ESPECIAL: C0: %f  C1: %f\n",temp1.v[0],temp1.v[1]);

                    ci0.v[ci0.x -1] = temp1.v[0];
                    ci1.v[ci1.x -1] = temp1.v[1];

                    //-------liberacion de memoria---------------
                    freev(temp1);   freev(tempX);   freev(tempHW);
                }


                //INTERPOLACION LINEAL DE VECTORES********
                temp1 = interp1(ti, ci0, t);
                copyRow(&C0,temp1,i);

                temp2 = interp1(ti, ci1, t);
                copyRow(&C1,temp2,i);

                temp3 = distancev(temp1,temp2);
                copyRow(&D,temp3,i);
                //***************************************

                //liberación de memoria
                freev(temp1);
                freev(temp2);
                freev(temp3);
                freev(ti);
                freev(hannW);
                freev(ci0);
                freev(ci1);

                i = i + nthreads;

            }
	    tFinal = omp_get_wtime();
	    double duracion = tFinal - tInicio;
	    printf("Termino el hilo %d en %f segundos\n",tid,duracion);
        //***********************************************************************************************************************
        } /* All threads join master thread and disband */
        //FIN PARTE 2
        //#############################################################################################################
        //#############################################################################################################
        //PARTE 3 DETERMINACION DEL UMBRAL DE CLASIFICACION Y POST-PROCESAMIENTO
	
	freev(xnp);
        //Eleccion de los mejores centroides para cada muestra en el tiempo, calculo de su punto medio y umbralizacion
        signal = classifySig(&x,&C0,&C1,&D);

        //Post-Procesamiento =)
        postProcesing(&x,&signal,BIAS);

        //Liberar toda la memoria utilizada
        freem(D);
        freem(C0);
        freem(C1);
        freev(t);
        freev(wsize);
        freev(wduration);

    }

	return signal;
}

//Weigthed 2 means
//RECIBE:   vector de la senal 'x' y vector de pesos 'w'
//REALIZA:  calcula los 2 centroides de las 2 clases obtenidas para la senal x con pesos w
//RETORNA:  devuelve un vector con los 2 centroides obtenidos
vec w2means(vec x, vec w){

    vec centroides = makev(2);  //clasificacion y centroides de las clases
    vec clasificacion = makev(x.x);  //vector de clasificacion
    vec clasifVieja = makev(x.x);    //vector de clasificacion dela iteracion anterior
    int maxIter = 100;      //numero maximo de iteraciones
    int i,j;                  //indiceS genericoS

    double sumc0, sumc1, sumw0, sumw1;

    centroides.v[0] = x.v[minv(x)]; //centroide de la clase 0 minima
    centroides.v[1] = x.v[maxv(x)]; //centroide de la clase 1 maxima

    /*printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n-*--*-*-*-*-*-*-*-*-*-*-*-***-*-\n\n");
    printf("        W2MEANS ALGORITHM        \n\n");
    printf("\n c0 = %f c1 = %f       \n\n",centroides.v[0],centroides.v[1]);
    printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n-*--*-*-*-*-*-*-*-*-*-*-*-***-*-\n\n");*/
    //system("PAUSE");

    //clasificaion inicial de x
    for (i = 0;i < x.x ;i++) {
        clasificacion.v[i] = (fabs((double)centroides.v[1]-x.v[i]) < fabs((double)x.v[i]-centroides.v[0])) ? 1:0;
    }

    i =0;
    clasifVieja.v[0] = -1;  //para que la primer vuelta sean distintos
    while (!equalv(clasifVieja,clasificacion) && i < maxIter){

       //printf("\n ''''''''''' iteracion w2means # %f",i);

       i++;
       freev(clasifVieja);
       clasifVieja = copyv(clasificacion);

        sumw0 = 0;
        sumw1 = 0;
        sumc0 = 0;
        sumc1 = 0;

        //Compute new weighted means
        for(j = 0;j < x.x ;j++){
            if (clasificacion.v[j]){ //clase 1
                sumw1 = (double) sumw1 + w.v[j];
                sumc1 = (double) sumc1 + x.v[j]*w.v[j] ;
            }
            else{                   //clase 0
                sumw0 = (double) sumw0 + w.v[j] ;
                sumc0 =  (double) sumc0 + x.v[j]*w.v[j];
            }
        }


        if (sumw1 == 0){sumw1 = 1;} //Verificacion para no dividir entre 0
        centroides.v[1] = (double) sumc1 / sumw1;
        centroides.v[0] = (double) sumc0 / sumw0;

        //NUEVA CLASIFICACION DE X
        for (j = 0;j < x.x ;j++) {
            clasificacion.v[j] = (fabs(x.v[j]-centroides.v[1]) < fabs(x.v[j]-centroides.v[0])) ? 1:0;
        }

    }

    /*printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n-*--*-*-*-*-*-*-*-*-*-*-*-***-*-\n\n");
    printf("        FIN DEL ALGORITMO W2MEANS ALGORITHM        \n\n");
    printf("\n c0 = %f c1 = %f       \n\n",centroides.v[0],centroides.v[1]);
    printf("-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n-*--*-*-*-*-*-*-*-*-*-*-*-***-*-\n\n");*/

    freev(clasificacion);
    freev(clasifVieja);

    return centroides;

}

//Algoritmo que calcula el vector de tamano de ventanas a usar
vec getWSvec (intervalo wdurationlim,double wincrement, double dt, int smax){

    int indice;
    double acumulador; //acumulador de los valores logaritmicos
    double limI; //limite en escala logaritmica
    double limS; //limite en escala logaritmica

    vec wsize;

   //VERIFICACION Y AJUSTE del limite inferior y superior del tamano en muestras de las ventanas*
    if (((double)wdurationlim.inicio/dt) >=_LOWERLIMIT_ && ((double)wdurationlim.inicio/dt) < smax) {
        limI = log2( (double)wdurationlim.inicio/dt );
    }
    else {//limite inferior invalido
        limI = log2(_LOWERLIMIT_);
        fprintf(stderr, "The Lower Window duration limit is invalid. Using default instead\n");
    }

    if (((double)wdurationlim.fin/dt) < smax && ((double)wdurationlim.fin/dt) > _LOWERLIMIT_){
        limS = log2( ((double)wdurationlim.fin/dt) );
    }
    else{//limite superior invalido
        limS = log2(smax -1);
        fprintf(stderr, "The Upper Window duration limit is invalid. Using default instead\n");
    }
    //*****************************************************************************************

    indice = 0;

    acumulador = limI;

    wsize =  zerov( ( (int)((limS-limI)/log2(wincrement)) ) + 1);

    //CASO PRIMERA VENTANA*************************************************************************************************
    //El caso de la primer ventana es especial, puede que al aproximarse a un multiplo de 2 se salga del intervalo definido
    if (((int)round(pow(2,acumulador))/2)*2 < pow(2,acumulador)) {//si la primer ventana no es valida
        acumulador = acumulador + log2(wincrement);
    }
    else{
        wsize.v[indice] = ((int)round(pow(2,acumulador))/2)*2; //solo multiplos de 2 en escala logaritmica
        indice = indice + 1;
        acumulador = acumulador + log2(wincrement);
    }
    //**********************************************************************************************************************

    //OBTENCION DE LAS VENTANAS*********************************************************************************************
    while (acumulador < limS) { //evaluacion en el plano logaritmico

        wsize.v[indice] = ((int)round(pow(2,acumulador))/2)*2; //solo multiplos de 2 en escala logaritmica

        if (indice < 1 || wsize.v[indice-1] != wsize.v[indice]) {
            indice = indice + 1;
        }
        acumulador = acumulador + log2(wincrement);
    }
    //********************************************************************************************************************

    //Resize al tamaño verdadero
    resizev(&wsize,indice);

    return wsize;
}

//Clasifica una señal de acuerdo a 2 matrices de centroides y una de distancias
vec classifySig(vec *x, matrix* C0,matrix* C1,matrix* D){

    vec classif = zerov(x->x);
    int i,j;
    int indOptimo;
    double umbral;

    //Para cada instante del tiempo
    for ( i = 0;i< D->y;i++){

        //Obtener los centroides mas separados
        indOptimo = 0;  //initial suposition

        for(j = 1;j < D->x; j++){
            if(D->m[j][i] > D->m[indOptimo][i]){
                indOptimo = j;
            }
        }

        //Umbralización
         umbral = (double)C1->m[indOptimo][i] - (double)(D->m[indOptimo][i]/2.);

         if (x->v[i] >= umbral){
            classif.v[i] = 1;
         }
    }
    return classif;
}

//funcion con valores por defecto
vec ts2means2(vec x){

    intervalo wdurationlim;
    wdurationlim.inicio = 4;
    wdurationlim.fin = x.x - 1;

    vec signal = ts2means(x,1,wdurationlim,sqrt(2),0.5,"max");

    return signal;
}

//funcion con valores por defecto, excepto el dt
vec ts2means3(vec x,double dt){

    intervalo wdurationlim;
    wdurationlim.inicio = 4*dt;
    wdurationlim.fin = (x.x - 1)*dt;
    
    if (wdurationlim.inicio >= wdurationlim.fin) fprintf(stderr,"Error: too short signal\n");
    vec signal = ts2means(x,dt,wdurationlim,sqrt(2),0.5,"max");

    return signal;
}

//funcion con valores por defecto, excepto el dt
void ts2means4(double * in,int sz,double dt,double * out){

    intervalo wdurationlim;
    wdurationlim.inicio = 4*dt;
    wdurationlim.fin = (sz - 1)*dt;
    
    //transform input double pointer into a vec type
    vec x = copyv2(in,sz);
    
    if (wdurationlim.inicio >= wdurationlim.fin) fprintf(stderr,"Error: too short signal\n");
      vec signal = ts2means(x,dt,wdurationlim,sqrt(2),0.5,"max");
    
    for (int i=0; i < signal.x; i++)
      out[i] = signal.v[i];
    
    //Could be done more efficiently  
    freev(signal);
    freev(x);

}


//Post-Procesado de la senal
void postProcesing(vec * x,vec * signal, char* BIAS){

    int j, i, iold, iChange, b;
    b = 0;

    if (strcmp(BIAS,"max") == 0){b=1;}
    intvec indicesCambios = makeiv(x->x + 1);

    for (i = 0;i < x->x +1; i++){
        indicesCambios.v[i] = -1;
    }

    //Determinar los indices de los cambios de clase
    iChange = 0;
    for (i = 1;i < x->x ;i++){
        if ((signal->v[i] - signal->v[i-1]) != 0){
            indicesCambios.v[iChange] = i;
            iChange++;
        }
    }

    i = 0;
    iold = 0;
    if (iChange != 0) { //Al menos debe haber un cambio -- BUG FIXED 2/5/2017 jcastro
	    //Ciclo de analisis de los cambios de clase----------------------------------------------------------------------
	    while (indicesCambios.v[i + 1] != -1) {

		//determine el tipo de cambio
		if((signal->v[indicesCambios.v[i]] - signal->v[indicesCambios.v[i] -1]) > 0)    {//de clase 0 a 1

		    if( minvv(getSubv(*x,iold, indicesCambios.v[i])) > maxvv(getSubv(*x,indicesCambios.v[i] +1, indicesCambios.v[i + 1])) ){
		        if (b) {
		            for(j = iold ;j <= indicesCambios.v[i];j++){
		                signal->v[j] = 1;
		            }
		        }
		        else{
		            for(j =  indicesCambios.v[i] +1; j <= indicesCambios.v[i + 1]; j++){
		                signal->v[j] = 0;
		            }
		        }
		    }

		}
		else {  //de clase 1 a 0

		    if( maxvv(getSubv(*x,iold, indicesCambios.v[i])) < minvv(getSubv(*x,indicesCambios.v[i] +1, indicesCambios.v[i +1])) ){
		        if (b) {
		            for(j =  indicesCambios.v[i] +1; j <= indicesCambios.v[i + 1]; j++){
		                signal->v[j] = 1;
		            }
		        }
		        else{
		            for(j = iold ;j <= indicesCambios.v[i];j++){
		                signal->v[j] = 0;
		            }
		        }
		    }
		}

		iold = indicesCambios.v[i];
		i++;
	    }
	    //fin del ciclo de los cambios de clases-------------------------------------------------------------------------

	   //ultimo caso-----------------------------------------------------------------------------------------------------
	    if(indicesCambios.v[i] != (x->x - 1)){

		if((signal->v[indicesCambios.v[i]] - signal->v[indicesCambios.v[i] -1]) > 0)    {//de clase 0 a 1

		    if( minvv(getSubv(*x,iold, indicesCambios.v[i])) > maxvv(getSubv(*x,indicesCambios.v[i] +1, (x->x - 1))) ){
		        if (b) {
		            for(j = iold ;j <= indicesCambios.v[i];j++){
		                signal->v[j] = 1;
		            }
		        }
		        else{
		            for(j =  indicesCambios.v[i] +1; j <= (x->x - 1); j++){
		                signal->v[j] = 0;
		            }
		        }
		    }

		}
		else {  //de clase 1 a 0

		    if( maxvv(getSubv(*x,iold, indicesCambios.v[i])) < minvv(getSubv(*x,indicesCambios.v[i] +1, (x->x - 1))) ){
		        if (b) {
		            for(j =  indicesCambios.v[i] +1; j <= (x->x - 1); j++){
		                signal->v[j] = 1;
		            }
		        }
		        else{
		            for(j = iold ;j <= indicesCambios.v[i]; j++){
		                signal->v[j] = 0;
		            }
		        }
		    }
		}

	    }
	    //-----------------------------------------------------------------------------------------------------------
    }
}
