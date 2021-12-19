//Voltage
//A = read("/home/naudin/Documents/article-2/AFD under Extreme Stimulation/Second Recordings/AFDnewDataSecondRecordingsNumberPointsDividedBy4.txt",-1,11);

//t=linspace(0,50,12501);

for i=[1:1:11]
//    plot2d(t,A(:,i),3)
end

//Steady-state current cost function
vecV=[-110:10:50]
Inf=[-68.6 -49.5 -18.2 -5.06 2.19 3.37 2.52 2.68 5.97 14.6 33.4 60.2 85 114 152 208 254]
//InfSD=[1 8.65 0.636 1.31 1.83 1.46 0.814 0.455 0.613 2.63 7.71 14.7 22.3 27.4 44.1 73.7 97.6]
InfSD=[1 8.65 0.636 1.31 1.83 0.01 0.01 0.455 0.613 2.63 7.71 14.7 22.3 27.4 44.1 73.7 97.6]

function y=WSS(pa)
    e=0;
    for i=1:length(vecV)
        tmp=0;
        tmp=(Inf(i)-(pa(1).*vecV(i).*vecV(i).*vecV(i) + pa(2).*vecV(i).*vecV(i) + pa(3).*vecV(i) + pa(4)))^2
        tmp=tmp/InfSD(i)
        e=e+tmp;
    end
    y=e/length(vecV) 
endfunction


function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=4;//nombre de paramètres à estimer
    pop=zeros(D,NP);//population : matrice de taille D*NP où NP est le nombre d'individus de ta population

    ////////////////////////////////////////////////////////////////////////
    //// Vecteurs de contraintes des paramètres : borne minimum/maximum ////
    ////////////////////////////////////////////////////////////////////////

    Xmin=[-200 -200 -200 -200];
    Xmax=[200  200  200  200];
    
    ///////////////////////////////////////////////////
    //// Initialisation aléatoire de ma population ////
    ///////////////////////////////////////////////////

    for j=1:NP
        for i=1:D
            pop(i,j)=Xmin(i)+(Xmax(i)-Xmin(i))*rand();
        end
    end
    
    //////////////////////////////////////////////////////////////
    //// Évaluation du meilleur individu après initialisation ////
    //////////////////////////////////////////////////////////////
    
    val=zeros(NP,1); // tableau avec le coût de chacun des individus, initialisé à 0
    
    // Evaluation de la fonction coût pour chacun des individus de la population
    for j=1:NP
        val(j)=WSS(pop(:,j))
    end
    
    ////////////////////////
    //// Étape suivante ////
    ////////////////////////
     
    iter=1; // nombre d'itération
    U=zeros(D,NP); // Vecteur intermédiaire perturbé (mutation + crossover)
    tempval=0;
    while iter<itermax
        for j=1:NP
            // ======= Construction de la matrice U = variation différentielle + crossover =======

            // ========= Tirage aléatoire de 3 entiers distincts r1, r2 et r3 et différents de j ========
            r1=j; r2=j; r3=j;//////////////////////////////////////
            while (r1==r2 | r1==r3 | r2==r3 | r1==j | r2==j | r3==j)
                r1=floor(1+NP*rand());
                r2=floor(1+NP*rand());
                r3=floor(1+NP*rand());
            end
            // ======== Variation différentielle =======
            V=pop(:,r1) + F*(pop(:,r2)-pop(:,r3));
            
            // ======== Contraintes ========
            for i=1:D
                if V(i)<=Xmin(i) then V(i)=Xmin(i);
                elseif V(i)>Xmax(i) then V(i)=Xmax(i);
                end
            end
            // ======== Crossover ========
            for i=1:D
                if rand()<CR then
                    U(i,j)=V(i);
                else
                    U(i,j)=pop(i,j);
                end
            end
        end // fin for j=1:NP
    
    // ======== Sélection ========
        for j=1:NP
            tempval=WSS(U(:,j));

            if tempval<=val(j) then
                pop(:,j) = U(:,j);
                val(j) = tempval;
            end
        end
        disp(iter) // c'est juste pour voir l'avancée de l'algorithme, à quelle itération on est...
        iter = iter + 1;
    end  //fin de la boucle while
    
    // Détermination de l'indice du meilleur individu
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    valBest=val(bestIndex);
    
    // Sauvegarde du meilleur individu
    bM = [];
    bM = pop(:,bestIndex);
    
    disp(val);
    disp(bM)
endfunction

//simulation(50,400,0.5,0.9)

for i=1:length(vecV)
    plot(vecV(i),Inf(i),'go')
    errbar(vecV(i),Inf(i),InfSD(i),InfSD(i))
end

function y=cubique(pa,V)
    y = pa(1).*V.*V.*V + pa(2).*V.*V + pa(3).*V + pa(4)
endfunction

pa=[   0.0003274
   0.0481816
   2.3119033
   38.98904]
   
//vecV=[-100:0.1:50]
//plot2d(vecV,cubique(pa,vecV),2)

for i=1:length(vecV)
    plot(vecV(i),cubique(pa,vecV(i)),'bx')
end

xlabel("$V$", "fontsize", 5)
ylabel("$I_{\infty}(V)$", "fontsize", 5)
a=gca();
a.font_size=4;
