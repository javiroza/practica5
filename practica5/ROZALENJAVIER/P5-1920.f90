! ---------------------------------- Pre-pràctica 5 ------------------------------------- !
! Autor: Javier Rozalén Sarmiento
! Grup: B1B
! Data: 05/11/2019
!
! Funcionalitat: es programen els mètodes de "sampleig" de densitats de probabilitat de canvi
! de variable i d'acceptació i rebuig, i es posen a prova amb algunes funcions.
!
! Comentaris: tot i que el meu grup és el B1B, tal i com vaig quedar amb en Bruno, entrego
! aquesta pràctica el mateix dia 05/11/2019 però a la sessió de les 15h-17h.

program pre_practica5
    implicit none
    double precision pi,M,a,b,fun,boxsize,valmitj,var,desvest
    double precision, allocatable :: xnums(:),xhis(:),vhis(:),errhis(:),xexpo(:)
    double precision, allocatable :: xhis2(:),vhis2(:),errhis2(:)
    integer ndat,i,ncaixes,ierr
    external fun
    common/cts/pi
    pi=acos(-1.d0)
    call srand(20034276)

    ! -------------------------------- Execici 1 --------------------------------------- !
    ! Assignem números i dimensions a les variables
    ndat=2000
    a=0.d0
    b=pi
    M=1.d0
    ncaixes=30
    allocate(xnums(ndat))
    allocate(xhis(ndat))
    allocate(vhis(ndat))
    allocate(errhis(ndat))

    ! Cridem les subrutines acceptrebuig i HISTOGRAMA per generar un arxiu de dades
    call acceptrebuig(ndat,xnums,a,b,M,fun)
    call HISTOGRAMA(ndat,xnums,a,b,ncaixes,xhis,vhis,errhis,boxsize,ierr)
    open(10,file="P5-1920-res.dat")
    write(10,*) "#xhis,vhis,tamany_caixa,barra_error"
    do i=1,ndat
        write(10,*) xhis(i),vhis(i),boxsize,errhis(i)
    enddo

    ! -------------------------------- Execici 2 --------------------------------------- !
    ndat=14000
    allocate(xexpo(ndat))
    allocate(xhis2(ndat))
    allocate(vhis2(ndat))
    allocate(errhis2(ndat))
    call subexpo(ndat,pi,xexpo)

    ! Càlcul del valor mitjà
    valmitj=0.d0
    do i=0,ndat
        valmitj=valmitj+xexpo(i)
    enddo
    valmitj=valmitj/dble(ndat)

    ! Càlcul de la variància
    var=0.d0
    do i=0,ndat
        var=var+(xexpo(i)-valmitj)**2.d0
    enddo
    var=var/dble(ndat)

    ! Càlcul de la desviació estàndard
    desvest=var**0.5d0

    print*,""
    print*,"Resultats exercici 2:"
    print*,valmitj,var,desvest

    ! Valor mitjà exacte: 1/pi = 0.31831
    ! Variància exacta:   1/pi² = 0.10132
    ! Desviació estàndard exacta: 1/pi = 0.31831

    ! Escriptura en fitxer dels valors que s'acaben de calcular
    write(10,*) ""
    write(10,*) "Valor mitjà calculat, Valor mitjà exacte",valmitj,"0.31831"
    write(10,*) "Variància calculada, variància exacta:",var,"0.10132"
    write(10,*) "Desviació estàndard calculada, desviació estàndard exacta:",desvest,"0.31831"
    write(10,*) ""
    write(10,*) ""

    ! Generació del segon histograma
    ncaixes=120
    a=0.d0
    b=3.d0

    open(11,file="aux.dat")
    call HISTOGRAMA(ndat,xexpo,a,b,ncaixes,xhis2,vhis2,errhis2,boxsize,ierr)
    do i=1,ndat
        write(11,*) xhis(i),vhis(i),boxsize,errhis(i)
    enddo

    close(10)
    close(11)
end program pre_practica5

! Subrutina acceptrebuig --> genera un vector amb ndat números distribuïts segons funci
subroutine acceptrebuig(ndat,xnums,a,b,M,funci)
    ! ndat --> nombre de números aleatoris que es vol generar (input)
    ! xnums --> vector que contindrà els ndat números aleatoris (output)
    ! a,b --> extrems de l'interval on està definida la nova variable aleatòria (input)
    ! M --> cota superior (input)
    ! funci --> densitat de probabilitat segons la qual està distribuida la nova var. aleat. (input)
    implicit none
    double precision xnums(ndat),a,b,M,funci,x1,x2,p,x
    double precision valmitj,var,desvest
    integer ndat,iseed,counter,i
    counter=0 ! Variable que porta el compte dels nombres aleatoris generats

    ! Generació dels nombres aleatoris
 1  x1=rand()
    x2=rand()
    x=(b-a)*x1+a
    p=M*x2
    if (funci(x).ge.p) then 
        xnums(counter)=x
        counter=counter+1
        if (counter.lt.ndat) then 
            goto 1
        endif
    else 
        goto 1
    endif

    ! Càlcul del valor mitjà
    valmitj=0.d0
    do i=0,ndat
        valmitj=valmitj+xnums(i)
    enddo
    valmitj=valmitj/dble(ndat)

    ! Càlcul de la variància
    var=0.d0
    do i=0,ndat
        var=var+(xnums(i)-valmitj)**2.d0
    enddo
    var=var/dble(ndat)

    ! "Càlcul" de la desviació estàndard
    desvest=var**0.5d0

    print*,""
    print*,"Resultats exercici 1:"
    print*,valmitj,var,desvest
    return
end subroutine acceptrebuig

! Funció fun --> distribució donada a l'exercici 1
double precision function fun(x)
    implicit none
    double precision pi,x
    common/cts/pi
    fun=(12.d0/(pi*(2.d0*pi**2.d0-3.d0)))*x**2.d0*(sin(x))**2.d0
    return
end function fun

! Subrutina subexpo --> genera un vector amb ndat números distribuïts exponencialment
subroutine subexpo(ndat,xlambda,xexpo)
    ! ndat --> nombre de números aleatoris (input)
    ! xlambda --> paràmetre de la distribució exponencial
    ! xexpo --> vector amb els nombres aleatoris exponencials
    implicit none
    integer ndat,i
    double precision xlambda,xexpo(ndat)

    do i=0,ndat
        xexpo(i)=(-1.d0/xlambda)*log(1-rand())
    enddo

    return
end subroutine subexpo

! Subrutina HISTOGRAMA --> agafada del campus virtual, genera un histograma
SUBROUTINE HISTOGRAMA(NDAT,XDATA,XA,XB,NBOX,XHIS,VHIS,ERRHIS,BOXSIZE,IERR)
    IMPLICIT NONE
    ! INPUT/OUTPUT VARIABLES
    INTEGER NDAT,NBOX
    DOUBLE PRECISION XDATA(NDAT),XA,XB
    DOUBLE PRECISION XHIS(NBOX),VHIS(NBOX),ERRHIS(NBOX)
    INTEGER IERR

    INTEGER I,IBOX,ICOUNT
    DOUBLE PRECISION BOXSIZE

    IF (XA.GE.XB) THEN 
        IERR=1
        RETURN
    ENDIF
    ! BOX SIZE
    BOXSIZE=(XB-XA)/NBOX

    ! COUNTS NUMBER OF POINTS WITHIN THE INTERVAL XA,XB
    ICOUNT=0

    ! SETS ALL TO ZERO
    DO I=1,NBOX
        VHIS(I)=0
        ERRHIS(I)=0
    ENDDO

    ! WE RUN THROUGH THE DATASET
    DO I=1,NDAT
    ! CHECKS IF DATA LIES WITHIN XA,XB
    IF (XDATA(I).GE.XA.AND.XDATA(I).LE.XB) THEN 
        IBOX=INT((XDATA(I)-XA)/BOXSIZE)+1
    ! PUTS XB INTO THE LAST BOX, IF NEEDED
        IF (IBOX.EQ.NBOX+1) IBOX=NBOX 

            VHIS(IBOX)=VHIS(IBOX)+1
            ICOUNT=ICOUNT+1
    ENDIF
    ENDDO

    IF (ICOUNT.EQ.0) THEN 
       IERR=2
       RETURN
    ENDIF

    IERR=0
    PRINT*,"ACCEPTED:",ICOUNT," OUT OF:",NDAT

    DO I=1,NBOX
    ! CENTRAL VALUE OF THE BAR
        XHIS(I)=XA+BOXSIZE/2.D0+(I-1)*BOXSIZE
    !  ERROBAR, STANDARD DEVIATION OF CORRESPONDING BINOMIAL
        ERRHIS(I)=SQRT(VHIS(I)/ICOUNT*(1.D0-VHIS(I)/ICOUNT))/BOXSIZE / SQRT(DBLE(ICOUNT))
    ! NORMALIZED VALUE OF THE BAR
    VHIS(I)=VHIS(I)/ICOUNT/BOXSIZE
    ENDDO
END SUBROUTINE HISTOGRAMA

