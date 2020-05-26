using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using FMatVecClr;

namespace EigvalClr
{
    class EigvalClrClass
    {
        static void Main(string[] args)
        {
            // Beispiel mit CLR arrays and statischer Methode
            double[] matrix = { 2.0, 30.0, 20.0, 4.0  };
            double[] eigValVec = LinearAlgebraDouble.eigval(2, matrix);
            foreach (double ev in eigValVec)
            {
                Console.WriteLine(ev);
            }

            // Beispiel mit Matrix Klasse und array return
            SquareMatrix sm = new SquareMatrix(2, matrix);
            eigValVec = sm.eigval();
            foreach (double ev in eigValVec)
            {
                Console.WriteLine(ev);
            }

            // Beispiel mit Matrix Klasse und ComplexVac return
            ComplexVec complexVec = sm.eigvalComplexVec();
            for(int i = 0; i < complexVec.Size; ++i)
            {
                Console.WriteLine(complexVec.e(i));
            }


        }
    }
}
