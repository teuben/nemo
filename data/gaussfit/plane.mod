parameter zero, gradx, grady;
data x,y;
observation z;

main()
{
    while(import())
        export(z - (zero + gradx*x + grady*y));
}


