#include "brains/SnapshotTestCase.hpp"
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SnapshotTestCase );

void SnapshotTestCase::setUp() {

}

void SnapshotTestCase::tearDown() {

}

void SnapshotTestCase::testMemoryLayout(){
    //test memory layout
    vector<Vector3d> v;
    const int num = 10000000;
    v.insert(v.end(), num, Vector3d(2.3, 0.84, 0.18));    

    double* sbegin = v[0].getArrayPointer();
    double* send = v[num-1].getArrayPointer() + 3;
    
    CPPUNIT_ASSERT_EQUAL(sizeof(Vector3d), sizeof(double) *3);    
    CPPUNIT_ASSERT_EQUAL((unsigned  int) sbegin,(unsigned  int) &v[0]);
    CPPUNIT_ASSERT_EQUAL((unsigned  int) send , (unsigned  int)&v[num]);

    //test memory access
    sbegin[12] = 32.01243;
    sbegin[13] = 1.023343;
    sbegin[14] = 82.025568;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(v[4][0], 32.01243, 0.00001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v[4][1], 1.023343, 0.00001);  
    CPPUNIT_ASSERT_DOUBLES_EQUAL(v[4][2], 82.025568, 0.00001);

    

}

void SnapshotTestCase::testConstructors(){
    Snapshot s;

    s.atomData.zAngle.push_back(1.0);
    
    double *p = Snapshot::getArrayPointer( s.atomData.zAngle);

}