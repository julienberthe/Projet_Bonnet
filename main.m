%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%              Projet Problemes Inverses                %
%                                                       %
%             Luc Laurent & Julien Berthe               %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

	%parametres materiaux
	donnee.mat.L=1000;			%longueur:
	donnee.mat.S=10^(-4);		%section:
	donnee.mat.E=220*10^9;		%module de young
	donnee.mat.rho=7000;		%masse volumique
    
    %Parametres de la methode EF
	donnee.nelem = 30;	%nombre d'elements
    
    %Approche inverse
    nbmode=2;   %nombre de modes propres observés à prendre en compte
    
    
    %mise en donnee
	for i=1:donnee.nelem
		donnee.Elem{i}.xinit=(i-1)*donnee.mat.L/donnee.nelem;
		donnee.Elem{i}.xfinal=i*donnee.mat.L/donnee.nelem;
		donnee.Elem{i}.dx=donnee.mat.L/donnee.nelem;
		donnee.Elem{i}.S=donnee.mat.S;
		donnee.Elem{i}.young=donnee.mat.E;
		donnee.Elem{i}.rho=donnee.mat.rho;
    end
    
    	donnee.dx =donnee.mat.L/donnee.nelem;
        donnee.x  =[0:donnee.dx:donnee.mat.L];
        
    %construction des differentes matrice du probleme EF
    disp('I	Construction des matrices');
    matrice=Construction_EF(donnee);
    
    %Prise en comtpe des CL (methode de substitution)
    disp('Ib   Conditionnement des matrices par subsitution')
    [matriceS,donneeS]=Substitution(matrice,donnee);
    
    disp('II   Résolution sur la base des modes propres')
    disp('  a	 Calcul des modes et valeurs propres');
    ModePropre=CalculModePropre(matriceS,donneeS);
    
    disp('  b    Normalisation des vecteurs propres')
    ModePropreNorme=NormalisationMode(ModePropre,matriceS);
    
    % Construction de la matrice R
    
    disp('III  Construction de la matrice R')
    matrice.Rw=ConstructionRw(ModePropreNorme,donnee,nbmode);
    matrice.Ru=ConstructionRu(ModePropreNorme,donnee,nbmode);
       
    % Calcul des valeurs observées    
    disp('IV  Calcul des valeurs observées')
    
    X=rand(size(ModePropreNorme.Vecteur,2),1);
    Vecteurp.orig=roundn(X,-3);
    clear X;
    ModePropreObs=Calculobs(Vecteurp.orig,ModePropreNorme,matriceR);
    
    % Probleme inverse pour retrouver deltap à partir des valeurs propres
    demarche='GCJ';     %%choix de l'approche: 'GCJ' gradient conjugué, 'QUAD' minimisation d'une fonction quadratique
        
    disp('V  Résolution du pb inverse');
    switch demarche;
        case 'GCJ';
            disp('>>>>> Méthode basée sur un gradient conjugué');
            tol=10^-9; %tolérance (critère d'arrêt du gradient conjugué)

            %%Données non bruitées
            disp('>>>>> Données non bruitées');
            delta.worig=ModePropreObs.Valeur-ModePropreNorme.Valeur';
            Vecteurp.recons=cgs(matriceR.Rw,delta.worig,tol);
            ModePropreComp.Valeur= ModePropreObs.Valeur - matriceR.Rw*Vecteurp.recons;
            Vecteurp_diff.orig=Vecteurp.recons-Vecteurp.orig;
            Vecteurp_relat.orig=zeros(donnee.nelem,1);
            for i=1:donnee.nelem
                Vecteurp_relat.orig(i)=Vecteurp_diff.orig(i)/Vecteurp.orig(i);
            end
            x=1:1:donnee.nelem;
            figure;
            bar(x,Vecteurp_diff.orig)
            bar(x,Vecteurp_relat.orig)

            % On bruite les données

            % Entre 0 et 1%
            disp('>>>>> Données bruitées (0% > 1%)');
            ModePropreObs.Bruit1.Valeur = ModePropreObs.Valeur*(1+10^-2*rand(1,1));
            delta.w1=ModePropreObs.Bruit1.Valeur-ModePropreNorme.Valeur';
            Vecteurp.recons_bruit1=cgs(matriceR.Rw,delta.w1,tol);
            Vecteurp_diff.bruit1=Vecteurp.recons_bruit1-Vecteurp.orig;
            Vecteurp_relat.bruit1=zeros(donnee.nelem,1);
            for i=1:donnee.nelem
                Vecteurp_relat.bruit1(i)=Vecteurp_diff.bruit1(i)/Vecteurp.orig(i);
            end

            % Entre 0 et 5%
            disp('>>>>> Données bruitées (0% > 5%)');
            ModePropreObs.Bruit2.Valeur = ModePropreObs.Valeur*(1+5*10-2*rand(1,1));
            delta.w2=ModePropreObs.Bruit2.Valeur-ModePropreNorme.Valeur';
            Vecteurp.recons_bruit2=cgs(matriceR.Rw,delta.w2,tol);
            Vecteurp_diff.bruit2=Vecteurp.recons_bruit2-Vecteurp.orig;
            Vecteurp_relat.bruit2=zeros(donnee.nelem,1);
            for i=1:donnee.nelem
                Vecteurp_relat.bruit2(i)=Vecteurp_diff.bruit2(i)/Vecteurp.orig(i);
            end

            % Entre 0 et 10%
            disp('>>>>> Données bruitées (0% > 10%)');
            ModePropreObs.Bruit3.Valeur = ModePropreObs.Valeur*(1+10^-1*rand(1,1));
            delta.w3=ModePropreObs.Bruit3.Valeur-ModePropreNorme.Valeur';
            Vecteurp.recons_bruit3=cgs(matriceR.Rw,delta.w3,tol);
            Vecteurp_diff.bruit3=Vecteurp.recons_bruit3-Vecteurp.orig;
            Vecteurp_relat.bruit3=zeros(donnee.nelem,1);
            for i=1:donnee.nelem
                Vecteurp_relat.bruit3(i)=Vecteurp_diff.bruit3(i)/Vecteurp.orig(i);
            end

            % Entre 0 et 20%
            disp('>>>>> Données bruitées (0% > 20%)');
            ModePropreObs.Bruit4.Valeur = ModePropreObs.Valeur*(1+0.2*rand(1,1));
            delta.w4=ModePropreObs.Bruit4.Valeur-ModePropreNorme.Valeur';
            Vecteurp.recons_bruit4=cgs(matriceR.Rw,delta.w4,tol);
            Vecteurp_diff.bruit4=Vecteurp.recons_bruit4-Vecteurp.orig;
            Vecteurp_relat.bruit4=zeros(donnee.nelem,1);
            for i=1:donnee.nelem
                Vecteurp_relat.bruit4(i)=Vecteurp_diff.bruit4(i)/Vecteurp.orig(i);
            end

            % Entre 0 et 50%
            disp('>>>>> Données bruitées (0% > 50%)');
            ModePropreObs.Bruit5.Valeur = ModePropreObs.Valeur*(1+0.5*rand(1,1));
            delta.w5=ModePropreObs.Bruit5.Valeur-ModePropreNorme.Valeur';
            Vecteurp.recons_bruit5=cgs(matriceR.Rw,delta.w5,tol);
            Vecteurp_diff.bruit5=Vecteurp.recons_bruit5-Vecteurp.orig;
            Vecteurp_relat.bruit5=zeros(donnee.nelem,1);
            for i=1:donnee.nelem
                Vecteurp_relat.bruit5(i)=Vecteurp_diff.bruit5(i)/Vecteurp.orig(i);
            end
            
        
        case 'QUAD';
             disp('>>>>> Méthode basée sur la minimisation d''une fonction quadratique');
            tol=10^-9; %tolérance 
            %%Données non bruitées
            disp('>>>>> Données non bruitées');
            Vecteurp.recons=min_quad();
            
    end
    %%Affichage
    Vect_comp=zeros(donnee.nelem,6);
    Vect_comp(:,1)=Vecteurp_diff.orig;
    Vect_comp(:,2)=Vecteurp_diff.bruit1;
    Vect_comp(:,3)=Vecteurp_diff.bruit2;
    Vect_comp(:,4)=Vecteurp_diff.bruit3;
    Vect_comp(:,5)=Vecteurp_diff.bruit4;
    Vect_comp(:,6)=Vecteurp_diff.bruit5;
    Vect_relat=zeros(donnee.nelem,6);
    Vect_relat(:,1)=Vecteurp_relat.orig;
    Vect_relat(:,2)=Vecteurp_relat.bruit1;
    Vect_relat(:,3)=Vecteurp_relat.bruit2;
    Vect_relat(:,4)=Vecteurp_relat.bruit3;
    Vect_relat(:,5)=Vecteurp_relat.bruit4;
    Vect_relat(:,6)=Vecteurp_relat.bruit5;
    
    
    figure;
    %plot(x,Vecteurp.recons,'b',x,Vecteurp.orig,'g',x,Vecteurp.recons_bruit1,'r');%,x,Vecteurp_recons_bruit1,'c',x,Vecteurp_recons_bruit2,'m',x,Vecteurp_recons_bruit3,'y',x,Vecteurp_recons_bruit4,'k');
    bar(x,Vect_comp,1.5)
    title('p=E/E0 reconstruit');
    figure;
    bar(x,Vect_relat,1.5)
    title('erreur relative sur p=E/E0 reconstruit');
    
    
    
    
    
    
    