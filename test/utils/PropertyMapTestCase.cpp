/*
 * Copyright (c) 2004-2020 The University of Notre Dame. All Rights Reserved.
 *
 * The University of Notre Dame grants you ("Licensee") a
 * non-exclusive, royalty free, license to use, modify and
 * redistribute this software in source and binary code form, provided
 * that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * This software is provided "AS IS," without a warranty of any
 * kind. All express or implied conditions, representations and
 * warranties, including any implied warranty of merchantability,
 * fitness for a particular purpose or non-infringement, are hereby
 * excluded.  The University of Notre Dame and its licensors shall not
 * be liable for any damages suffered by licensee as a result of
 * using, modifying or distributing the software or its
 * derivatives. In no event will the University of Notre Dame or its
 * licensors be liable for any lost revenue, profit or data, or for
 * direct, indirect, special, consequential, incidental or punitive
 * damages, however caused and regardless of the theory of liability,
 * arising out of the use of or inability to use software, even if the
 * University of Notre Dame has been advised of the possibility of
 * such damages.
 *
 * SUPPORT OPEN SCIENCE!  If you use OpenMD or its source code in your
 * research, please cite the appropriate papers when you publish your
 * work.  Good starting points are:
 *
 * [1] Meineke, et al., J. Comp. Chem. 26, 252-271 (2005).
 * [2] Fennell & Gezelter, J. Chem. Phys. 124, 234104 (2006).
 * [3] Sun, Lin & Gezelter, J. Chem. Phys. 128, 234107 (2008).
 * [4] Vardeman, Stocker & Gezelter, J. Chem. Theory Comput. 7, 834 (2011).
 * [5] Kuang & Gezelter, Mol. Phys., 110, 691-701 (2012).
 * [6] Lamichhane, Gezelter & Newman, J. Chem. Phys. 141, 134109 (2014).
 * [7] Lamichhane, Newman & Gezelter, J. Chem. Phys. 141, 134110 (2014).
 * [8] Bhattarai, Newman & Gezelter, Phys. Rev. B 99, 094106 (2019).
 */

#include "utils/PropertyMapTestCase.hpp"
#include <iostream>
#include <algorithm>
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( PropertyMapTestCase );


void PropertyMapTestCase::testPropertyMap(){
    PropertyMap props;

    //test addProperty
    BoolGenericData* b0 = new BoolGenericData("BoolData");
    b0->setData(false);
    props.addProperty(b0);
    CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b0);
    
    BoolGenericData* b1 = new BoolGenericData("BoolData");
    b1->setData(true);
    props.addProperty(b1);
    CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b1);


    IntGenericData* i1 = new IntGenericData("IntData");
    i1->setData(89);
    props.addProperty(i1);
    
    FloatGenericData* f1 = new FloatGenericData("FloatData");
    f1->setData(49.328);
    props.addProperty(f1);

    DoubleGenericData* d1 = new DoubleGenericData("DoubleData");
    d1->setData(95.1933432);
    props.addProperty(d1);

    StringGenericData* s1 = new StringGenericData("StringData");
    s1->setData("Hello");
    props.addProperty(s1);


    IntVectorGenericData* iv1 = new IntVectorGenericData("IntVector");
    iv1->push_back(2);
    iv1->push_back(1);
    iv1->push_back(324);
    props.addProperty(iv1);

    //test getPropertyNames
    std::vector<std::string> propNames = props.getPropertyNames();


    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") != propNames.end());

    //test getProperties    
    std::vector<GenericData*> propPointers = props.getProperties();    
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), b1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), i1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), f1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), d1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), s1) != propPointers.end());
    CPPUNIT_ASSERT(std::find(propPointers.begin(), propPointers.end(), iv1) != propPointers.end());

    //test getPropertyByName
    CPPUNIT_ASSERT(props.getPropertyByName("BoolData") == b1);
    CPPUNIT_ASSERT(props.getPropertyByName("IntData") == i1);
    CPPUNIT_ASSERT(props.getPropertyByName("FloatData") == f1);
    CPPUNIT_ASSERT(props.getPropertyByName("DoubleData") == d1);
    CPPUNIT_ASSERT(props.getPropertyByName("StringData") == s1);
    CPPUNIT_ASSERT(props.getPropertyByName("IntVector") == iv1);

    CPPUNIT_ASSERT(b1->getData() == true);
    CPPUNIT_ASSERT(i1->getData() == 89);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f1->getData(), 49.328, 0.000001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(d1->getData(), 95.1933432, 0.000001);   
    CPPUNIT_ASSERT(s1->getData() == "Hello");
    CPPUNIT_ASSERT_EQUAL((*iv1)[0], 2);
    CPPUNIT_ASSERT_EQUAL((*iv1)[1], 1);
    CPPUNIT_ASSERT_EQUAL((*iv1)[2], 324);        

    //test removeProperty
    props.removeProperty("DoubleData");
    props.removeProperty("FloatData");
    props.removeProperty("IntVector");
    propNames = props.getPropertyNames();
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") != propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") == propNames.end());

    //test clearProperties
    props.clearProperties();
    propNames = props.getPropertyNames();
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "BoolData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "FloatData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "DoubleData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "StringData") == propNames.end());
    CPPUNIT_ASSERT(std::find(propNames.begin(), propNames.end(), "IntVector") == propNames.end());    

}
