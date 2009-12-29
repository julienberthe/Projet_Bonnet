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
	donnee.nelem = 16;	%nombre d'elements
    
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
    
    disp('II   R√©solution sur la base des modes propres')
    disp('  a	 Calcul des modes et valeurs propres');
    ModePropre=CalculModePropre(matriceS,donneeS);
    
    disp('  b    Normalisation des vecteurs propres')
    ModePropreNorme=NormalisationMode(ModePropre,matriceS);
    
    % Construction de la matrice R
    
    disp('III  Construction de la matrice R')
    matriceR1=ConstructionRw(ModePropreNorme,matrice,donnee);
    matriceR=ConstructionRu(ModePropreNorme,matrice,donnee,matriceR1);
    matriceR.Rw=matriceR1.Rw;
    clear matriceR1;
    
    % Calcul des valeurs observ√©es
    
    disp('IV  Calcul des valeurs observ√©es')
    X=rand(size(ModePropreNorme.Vecteur,2),1);
    Vecteurp=roundn(X,-3);
    clear X;
    ModePropreObs=Calculobs(Vecteurp,ModePropreNorme,matriceR);
    
    % Probleme inverse pour retrouver deltap ‡ partir des valeurs propres
    
    deltaw=ModePropreObs.Valeur-ModePropreNorme.Valeur';
    Vecteurp_recons=cgs(matriceR.Rw,deltaw);
    ModePropreComp.Valeur= ModePropreObs.Valeur - matriceR.Rw*Vecteurp_recons;
    x=1:1:size(ModePropreNorme.Valeur,2);
    
    % On bruite les donnÈes
    
    % Entre 0 et 1%
   
    ModePropreObsBruit1.Valeur = ModePropreObs.Valeur + (ModePropreObs.Valeur/100)*rand(1,1);
    deltaw1=ModePropreObsBruit1.Valeur-ModePropreNorme.Valeur';
    Vecteurp_recons_bruit1=cgs(matriceR.Rw,deltaw1);
    
    
    % Entre 0 et 5%
    
    ModePropreObsBruit2.Valeur = ModePropreObs.Valeur + (ModePropreObs.Valeur/20)*rand(1,1);
    deltaw2=ModePropreObsBruit2.Valeur-ModePropreNorme.Valeur';
    Vecteurp_recons_bruit2=cgs(matriceR.Rw,deltaw2);
    
    % Entre 0 et 10%
   
    ModePropreObsBruit3.Valeur = ModePropreObs.Valeur + (ModePropreObs.Valeur/10)*rand(1,1);
    deltaw3=ModePropreObsBruit3.Valeur-ModePropreNorme.Valeur';
    Vecteurp_recons_bruit3=cgs(matriceR.Rw,deltaw3);
    
    % Entre 0 et 20%
   
    ModePropreObsBruit4.Valeur = ModePropreObs.Valeur + (ModePropreObs.Valeur/5)*rand(1,1);
    deltaw4=ModePropreObsBruit4.Valeur-ModePropreNorme.Valeur';
    Vecteurp_recons_bruit4=cgs(matriceR.Rw,deltaw4);
    
    % Entre 0 et 100%
    
    ModePropreObsBruit5.Valeur = ModePropreObs.Valeur + (ModePropreObs.Valeur/1)*rand(1,1);
    deltaw5=ModePropreObsBruit5.Valeur-ModePropreNorme.Valeur';
    Vecteurp_recons_bruit5=cgs(matriceR.Rw,deltaw5);
    
    plot(x,Vecteurp_recons,'b',x,Vecteurp,'g',x,Vecteurp_recons_bruit1,'r');%,x,Vecteurp_recons_bruit1,'c',x,Vecteurp_recons_bruit2,'m',x,Vecteurp_recons_bruit3,'y',x,Vecteurp_recons_bruit4,'k');
    
    
    
    
    
    
    
    