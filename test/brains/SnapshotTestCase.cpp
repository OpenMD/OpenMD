#include "brains/SnapshotTestCase.hpp"
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SnapshotTestCase );

void SnapshotTestCase::setUp() {
}

void SnapshotTestCase::tearDown() {
}
void SnapshotTestCase::testConstructors(){
    Snapshot s;

    s.pos.push_back(Vector3d(0, 1, 2));
    s.zAngle.push_back(1.0);
    
    double *p = Snapshot::getArrayPointer(s.zAngle);

    p = Snapshot::getArrayPointer(s.pos);
}
