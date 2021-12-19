//Steady-state current
vecV=[-120:10:50];
Inf=[-13.1 -10.4 -7.92 -5.89 -4.11 -2.69 -1.02 0.0211 1.17 3.1 7.32 14.2 22.4 31.5 43.2 54.5 69.5 82.4];
InfSD=[2.88 2.55 1.47 1.31 1.04 0.809 0.7 0.658 0.638 0.889 1.94 3.5 5.36 7.63 10.6 13.3 16 17.9]

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
    y = pa(1).*vecV(i).*vecV(i).*vecV(i) + pa(2).*vecV(i).*vecV(i) + pa(3).*vecV(i) + pa(4)
endfunction

pa=[   0.0000438
   0.0093345
   0.7727631
   20.380413]

for i=1:length(vecV)
    plot(vecV(i),cubique(pa,vecV(i)),'bx')
end

xlabel("$V$", "fontsize", 5)
ylabel("$I_{\infty}(V)$", "fontsize", 5)
a=gca();
set(gca(),"data_bounds",matrix([-122,60,-20,120],2,-1));
a.font_size=4;

