using System;
using System.Runtime.InteropServices;

namespace eigval_net
{
    class EigVal
    {
        // wrapper dll must be in PATH!
        [DllImport(@"wrapper.dll", EntryPoint = "eigval", CallingConvention = CallingConvention.StdCall)]
        public static extern void eigval(UInt64 s, double[] squareMatrix, double[] eigValVec);

        static void Main(string[] args)
        {
            UInt64 s = 2;
            double[] squareMatrix = { 2.0, 30.0, 20.0, 4.0 };
            double[] eigValVec = new double[s];
            eigval(s, squareMatrix, eigValVec);
            foreach(double ev in eigValVec)
            {
                Console.WriteLine(ev);
            }

        }
    }
}
