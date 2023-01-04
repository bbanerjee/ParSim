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
  ps->output("-");

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
  out_ps->output("-");

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

  // Create manager
  try {
    MaterialManager manager;

    // Register material
    manager.registerMaterial("material 0", mat0);
    manager.registerMaterial("material 1", mat1);
    manager.registerMaterial("material 2", mat2, 7);
    manager.registerEmptyMaterial(mat3);

    manager.finalizeMaterials();
  } catch (const ProblemSetupException& e) {
    std::cout << e.message() << "\n";
  } catch (const std::exception& e) {
    std::cout << e.what() << "\n";
  }
}