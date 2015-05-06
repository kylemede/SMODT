#!/usr/bin/env python
import unittest
import math
import hello
import orbCalc
import numpy as np

class TestSwig(unittest.TestCase):
    """
    Test linking with swig. Use simple cases to mimic the functions 
    we will encounter in this package. If this test fail, then the 
    package is not portable on the local machine. 
    """
    def testInt(self):
        #pass a int from python as an argument of a c++ function, return an int
        self.assertEqual(hello.fact(4),24)
    def testFloat(self):
        #pass a float from python as an argument of a c++ function, return an
        #float.
        self.assertAlmostEqual(hello.csin(0.5*math.pi),math.sin(0.5*math.pi))
    def testStr(self):
        #pass a string from python as an argument of a c++ function.
        hello.speak("helloworld")
        self.assertAlmostEqual(1,1)
    def testbool(self):
        #pass a bool from python as an argument of a c++ function.
        a = hello.selection(True)
        self.assertEqual(a,1)
        b = hello.selection(False)
        self.assertEqual(b,0)
    def testNparray1d_scalar(self):
        #pass an 1-d numpy array from python as an argument of a c++ function.
        #the value of the array is not changed in the c++ function.
        a = np.arange(50)
        self.assertAlmostEqual(hello.sum_scalar(a,1),np.sum(a)+50)
    def testNparray1d_1d(self):
        #pass two 1-d numpy array from python as an argument of a c++ function.
        #the value of the array is not changed in the c++ function.
        a = np.arange(50)
        self.assertAlmostEqual(hello.sum_1d(a,a),np.sum(a)*2)
    def testNparray2d(self):
        #pass a 2-d numpy array from python as an argument of a c++ function.
        #the value of the array is not changed in the c++ function.
        a = np.array([[1,2],[3,4]])
        self.assertAlmostEqual(hello.sum2d(a),10)	

    # tests for initial version of orbCalc.cpp
    def testInt(self):
        #pass a int from python as an argument of a c++ function, return an int
        self.assertEqual(orbCalc.fact(4),24)
    def testFloat(self):
        #pass a float from python as an argument of a c++ function, return an
        #float.
        self.assertAlmostEqual(orbCalc.csin(0.5*math.pi),math.sin(0.5*math.pi))
    def testStr(self):
        #pass a string from python as an argument of a c++ function.
        orbCalc.speak("orbCalcworld")
        self.assertAlmostEqual(1,1)
    def testbool(self):
        #pass a bool from python as an argument of a c++ function.
        a = orbCalc.selection(True)
        self.assertEqual(a,1)
        b = orbCalc.selection(False)
        self.assertEqual(b,0)
    def testNparray1d_scalar(self):
        #pass an 1-d numpy array from python as an argument of a c++ function.
        #the value of the array is not changed in the c++ function.
        a = np.arange(50)
        self.assertAlmostEqual(orbCalc.sum_scalar(a,1),np.sum(a)+50)
    def testNparray1d_1d(self):
        #pass two 1-d numpy array from python as an argument of a c++ function.
        #the value of the array is not changed in the c++ function.
        a = np.arange(50)
        self.assertAlmostEqual(orbCalc.sum_1d(a,a),np.sum(a)*2)
    def testNparray2d(self):
        #pass a 2-d numpy array from python as an argument of a c++ function.
        #the value of the array is not changed in the c++ function.
        a = np.array([[1,2],[3,4]])
        self.assertAlmostEqual(orbCalc.sum2d(a),10)		
if __name__ == "__main__":
    unittest.main()

