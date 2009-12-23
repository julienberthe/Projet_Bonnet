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
    
    disp('II   Résolution sur la base des modes propres')
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
    
    % Calcul des valeurs observées
    
    disp('IV  Calcul des valeurs observées')
    X=rand(size(ModePropreNorme.Vecteur,2),1);
    Vecteurp=roundn(X,-3);
    clear X;
    ModePropreObs=Calculobs(Vecteurp,ModePropreNorme,matriceR);
    
    
    
    