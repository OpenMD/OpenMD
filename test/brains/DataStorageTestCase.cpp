/*
 * Copyright (c) 2004-present, The University of Notre Dame. All rights
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the following paper when you publish your work:
 *
 * [1] Drisko et al., J. Open Source Softw. 9, 7004 (2024).
 *
 * Good starting points for code and simulation methodology are:
 *
 * [2] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [3] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [4] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [5] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [6] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [7] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 * [9] Drisko & Gezelter, J. Chem. Theory Comput. 20, 4986-4997 (2024).
 */

#include "brains/DataStorageTestCase.hpp"
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(DataStorageTestCase);

void DataStorageTestCase::testDataStorage() {
  // test constructor
  DataStorage ds1(10020);
  CPPUNIT_ASSERT_EQUAL(ds1.getStorageLayout(), 0x11111111);
  CPPUNIT_ASSERT_EQUAL(ds1.getSize(), 10020);

  // test resize
  ds1.resize(56);
  CPPUNIT_ASSERT_EQUAL(ds1.getSize(), 56);

  // test reserve

  DataStorage ds2(10000, DataStorage::dslForce | DataStorage::dslAmat |
                             DataStorage::dslZAngle);

  CPPUNIT_ASSERT(!(ds2.getStorageLayout() & DataStorage::dslPosition));
  CPPUNIT_ASSERT(!(ds2.getStorageLayout() & DataStorage::dslVelocity));
  CPPUNIT_ASSERT(ds2.getStorageLayout() & DataStorage::dslAmat);
  CPPUNIT_ASSERT(!(ds2.getStorageLayout() & DataStorage::dslAngularMomentum));
  CPPUNIT_ASSERT(!(ds2.getStorageLayout() & DataStorage::dslUnitVector));
  CPPUNIT_ASSERT(ds2.getStorageLayout() & DataStorage::dslZAngle);
  CPPUNIT_ASSERT(ds2.getStorageLayout() & DataStorage::dslForce);
  CPPUNIT_ASSERT(!(ds2.getStorageLayout() & DataStorage::dslTorque));

  CPPUNIT_ASSERT(ds2.position.size() == 0);
  CPPUNIT_ASSERT(ds2.velocity.size() == 0);
  CPPUNIT_ASSERT(ds2.aMat.size() == 10000);
  CPPUNIT_ASSERT(ds2.angularMomentum.size() == 0);
  CPPUNIT_ASSERT(ds2.unitVector.size() == 0);
  CPPUNIT_ASSERT(ds2.zAngle.size() == 10000);
  CPPUNIT_ASSERT(ds2.force.size() == 10000);
  CPPUNIT_ASSERT(ds2.torque.size() == 0);

  // test getArrayPointer

  double* pamat   = ds2.getArrayPointer(DataStorage::dslAmat);
  double* pzangle = ds2.getArrayPointer(DataStorage::dslZAngle);
  double* pforce  = ds2.getArrayPointer(DataStorage::dslForce);

  pamat[9] = 2.0993;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[1](0, 0), 2.0993,
                               OpenMD::NumericConstant::epsilon);

  pzangle[9743] = 8464353.129574;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.zAngle[9743], 8464353.129574,
                               OpenMD::NumericConstant::epsilon);

  pforce[12] = 345546.3413263;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.force[4][0], 345546.3413263,
                               OpenMD::NumericConstant::epsilon);

  // test copy
  ds2.zAngle[0] = 7.0;
  ds2.zAngle[1] = 3.5;

  ds2.force[0][0] = 6;
  ds2.force[0][1] = 7.0;
  ds2.force[0][2] = 8.0;

  ds2.force[1][0] = 9.0;
  ds2.force[1][1] = 10.0;
  ds2.force[1][2] = 11.0;

  ds2.aMat[0](0, 0) = 1.0;
  ds2.aMat[0](0, 1) = 2.0;
  ds2.aMat[0](0, 2) = 3.0;
  ds2.aMat[0](1, 0) = 4.0;
  ds2.aMat[0](1, 1) = 5.0;
  ds2.aMat[0](1, 2) = 6.0;
  ds2.aMat[0](2, 0) = 7.0;
  ds2.aMat[0](2, 1) = 8.0;
  ds2.aMat[0](2, 2) = 9.0;

  ds2.aMat[1](0, 0) = 11.0;
  ds2.aMat[1](0, 1) = 12.0;
  ds2.aMat[1](0, 2) = 13.0;
  ds2.aMat[1](1, 0) = 14.0;
  ds2.aMat[1](1, 1) = 15.0;
  ds2.aMat[1](1, 2) = 16.0;
  ds2.aMat[1](2, 0) = 17.0;
  ds2.aMat[1](2, 1) = 18.0;
  ds2.aMat[1](2, 2) = 19.0;

  ds2.copy(0, 2, 10);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.zAngle[10], 7.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.zAngle[11], 3.5,
                               OpenMD::NumericConstant::epsilon);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.force[10][0], 6.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.force[10][1], 7.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.force[10][2], 8.0,
                               OpenMD::NumericConstant::epsilon);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.force[11][0], 9.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.force[11][1], 10.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.force[11][2], 11.0,
                               OpenMD::NumericConstant::epsilon);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](0, 0), 1.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](0, 1), 2.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](0, 2), 3.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](1, 0), 4.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](1, 1), 5.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](1, 2), 6.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](2, 0), 7.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](2, 1), 8.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[10](2, 2), 9.0,
                               OpenMD::NumericConstant::epsilon);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](0, 0), 11.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](0, 1), 12.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](0, 2), 13.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](1, 0), 14.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](1, 1), 15.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](1, 2), 16.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](2, 0), 17.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](2, 1), 18.0,
                               OpenMD::NumericConstant::epsilon);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ds2.aMat[11](2, 2), 19.0,
                               OpenMD::NumericConstant::epsilon);
}
