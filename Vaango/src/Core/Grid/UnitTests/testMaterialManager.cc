#include <Core/Exceptions/ProblemSetupException.h>
#include <Core/Grid/Material.h>
#include <Core/Grid/EmptyMaterial.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/ProblemSpec/ProblemSpec.h>
#include <gtest/gtest.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

#include <iostream>

using namespace Uintah;

TEST(MaterialTest, construction) {
  // Default
  Material mat0;

  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "material");
  xmlNewProp(rootNode, BAD_CAST "name", BAD_CAST "test material");
  xmlDocSetRootElement(doc, rootNode);

  // Print the document to stdout
  // xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }
  //ps->output("-");

  // With ProblemSpec
  Material mat(ps);

  // Check names
  EXPECT_EQ(mat.hasName(), true);

  // Get name
  auto name = mat.getName();
  EXPECT_EQ(name, "test material");

  // Set the matID
  mat.setDWIndex(5);
  EXPECT_ANY_THROW({ mat.setDWIndex(6); });

  // Get the matID
  int id = mat.getDWIndex();
  EXPECT_EQ(id, 5);

  // Output
  ProblemSpecP out_ps = mat.outputProblemSpec(ps);
  //out_ps->output("-");

  // Get the material subset
  auto mat_ptr = mat.thisMaterial();
  auto data    = mat_ptr->getVector();
  for (auto& idx : data) {
    EXPECT_EQ(idx, 5);
  }
}

TEST(MaterialManagerTest, construction) {
  // Basic constructor
  EXPECT_NO_THROW({
    MaterialManager manager;
    EXPECT_THROW({ MaterialManager manager1; }, ProblemSetupException);
  });
}

TEST(MaterialManagerTest, others) {
  // Create a new document
  xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode = xmlNewNode(nullptr, BAD_CAST "material");
  xmlNewProp(rootNode, BAD_CAST "name", BAD_CAST "test material");
  xmlDocSetRootElement(doc, rootNode);

  // Print the document to stdout
  // xmlSaveFormatFileEnc("-", doc, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps = scinew ProblemSpec(xmlDocGetRootElement(doc), false);
  if (!ps) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  // Create a material With ProblemSpec
  std::shared_ptr<Material> mat0 = std::make_shared<Material>(ps);
  std::shared_ptr<Material> mat1 = std::make_shared<Material>(ps);
  std::shared_ptr<Material> mat2 = std::make_shared<Material>(ps);
  std::shared_ptr<EmptyMaterial> mat3 = std::make_shared<EmptyMaterial>();
  std::shared_ptr<Material> mat4 = std::make_shared<Material>(ps);

  // Create manager
  try {
    MaterialManager manager;

    // Register material
    manager.registerMaterial("material 0", mat0);
    manager.registerMaterial("material 1", mat1);
    manager.registerMaterial("material 2", mat2, 7);
    manager.registerEmptyMaterial(mat3);
    manager.registerMaterial("material 1", mat4);

    manager.finalizeMaterials();

    // Get materials
    [[maybe_unused]] auto matset = manager.allMaterials();
    [[maybe_unused]] auto matset1 = manager.allMaterials("material 1");

    // Get num materials
    auto num = manager.getNumMaterials();
    auto num1 = manager.getNumMaterials("material 1");
    EXPECT_EQ(num, 5);
    EXPECT_EQ(num1, 2);
    //std::cout << "num materials: " << num << " matset: " << *matset << "\n";
    //std::cout << "num1 materials: " << num1 << " matset1: " << *matset1 << "\n";

    // Get materials
    auto mat_get00 = manager.getMaterial(15);
    EXPECT_EQ(mat_get00, nullptr);
    auto mat_get01 = manager.getMaterial("material 1", 7);
    EXPECT_EQ(mat_get01, nullptr);

    auto mat_get = manager.getMaterial(1);
    if (mat_get) {
      auto mat_get_ptr = mat_get->thisMaterial();
      auto data_get    = mat_get_ptr->getVector();
      for (auto& idx : data_get) {
        //std::cout << "mat_get: idx = " << idx << "\n";
        EXPECT_EQ(idx, 1);
      }
    }

    auto mat_get1 = manager.getMaterial("material 1", 1);
    if (mat_get1) {
      auto mat_get_ptr = mat_get1->thisMaterial();
      auto data_get    = mat_get_ptr->getVector();
      for (auto& idx : data_get) {
        //std::cout << "mat_get1: idx = " << idx << "\n";
        EXPECT_EQ(idx, 4);
      }
    }

    /* Now private*/
    //auto mat_get02 = manager.getMaterialByName("test material");

    // Get all materials in one
    [[maybe_unused]] auto allmat = manager.getAllInOneMaterial();
    //std::cout << "All = " << *allmat << "\n";

    // Get original materials
    [[maybe_unused]] auto orig = manager.originalAllMaterials();
    //std::cout << "Original = " << *orig << "\n";

  } catch (const ProblemSetupException& e) {
    std::cout << e.message() << "\n";
  } catch (const std::exception& e) {
    std::cout << e.what() << "\n";
  }

}

TEST(MaterialManagerTest, forICE) {

  // Create a new document
  xmlDocPtr doc1 = xmlNewDoc(BAD_CAST "1.0");

  // Create root node
  xmlNodePtr rootNode1 = xmlNewNode(nullptr, BAD_CAST "test");
  xmlDocSetRootElement(doc1, rootNode1);
  xmlNewChild(rootNode1, nullptr, BAD_CAST "from_material",
              BAD_CAST "material from");

  // Print the document to stdout
  //xmlSaveFormatFileEnc("-", doc1, "ISO-8859-1", 1);

  // Create a ProblemSpec
  ProblemSpecP ps1 = scinew ProblemSpec(xmlDocGetRootElement(doc1), false);
  if (!ps1) {
    std::cout << "**Error** Could not create ProblemSpec." << std::endl;
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
    exit(-1);
  }

  MaterialManager manager;
  std::shared_ptr<Material> mat0 = std::make_shared<Material>(ps1);
  manager.registerMaterial("from_material", mat0);
  manager.finalizeMaterials();

  ProblemSpecP child = ps1->findBlock("test");
  auto mat_get02 = manager.parseAndLookupMaterial(child, "from_material");
  if (mat_get02) {
    auto mat_get_ptr = mat_get02->thisMaterial();
    auto data_get    = mat_get_ptr->getVector();
    for (auto& idx : data_get) {
      //std::cout << "mat_get2: idx = " << idx << "\n";
      EXPECT_EQ(idx, 0);
    }
  }
}
