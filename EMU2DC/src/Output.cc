#include <Output.h>
#include <Exception.h>
#include <NodePArray.h>
#include <Body.h>
#include <Node.h>

#include <Core/ProblemSpec/ProblemSpec.h>

#include <unistd.h>
#include <fstream>

using namespace Emu2DC;

Output::Output()
{
}

Output::Output(const Uintah::ProblemSpecP& ps)
{
  initialize(ps);
}

Output::~Output()
{
}

void
Output::initialize(const Uintah::ProblemSpecP& ps)
{
  // Get output files/interval
  Uintah::ProblemSpecP io_ps = ps->findBlock("Output");
  if (!io_ps) {
    throw Exception("**ERROR** <Output> tag not found", __FILE__, __LINE__); 
  } 

  io_ps->require("output_file", d_output_file_name);
  io_ps->require("output_iteration_interval", d_output_iter_interval);
  char buffer[2000];
  char * str = getcwd( buffer, 2000 );
  if (str == NULL) {
    throw Exception("**ERROR** Directory not returned by getcwd()", __FILE__, __LINE__); 
  } else {
    d_output_folder_name = std::string(buffer);
  }
  d_output_file_count = 0;
}

void
Output::write(const Time& time, const BodySPArray& bodyList) 
{
  // Write the output to individual files
  std::string output_file_name = d_output_file_name;
  std::string current_output_file_name = d_output_file_name + std::to_string(d_output_file_count) + ".tec";
  //write(current_output_file_name, fmt='(A,I5.5,A)') trim(output_file_name), output_file_count, '.tec'

  std::ofstream output_file(current_output_file_name);

  // Loop through bodies
  for (auto body_iter = bodyList.begin(); body_iter != bodyList.end(); ++body_iter) {

    const NodePArray& node_list = (*body_iter)->nodes();

    // count valid nodes
    int valid_node_count = 0;
    for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {
      if ((*node_iter)->omit()) continue;  // skip this node
      valid_node_count++;
    }

    // write the headers
    output_file << "TITLE=\"simulation results\" " << std::endl;
    output_file << "VARIABLES=\"X\",\"Y\",\"DX\",\"DY\",\"VX\",\"VY\",\"DAM\",\"W\" " << std::endl;
    output_file << "ZONE I=" << valid_node_count << " SOLUTIONTIME=" << time.currentTime() 
                << " F=POINT" << std::endl;

    for (auto node_iter = node_list.begin(); node_iter != node_list.end(); ++node_iter) {
      NodeP cur_node = *node_iter;
      if (cur_node->omit()) continue;  // skip this node
      double xdisp = cur_node->displacement()[0];
      double ydisp = cur_node->displacement()[1];
      double cur_x_pos = cur_node->position()[0] + xdisp;
      double cur_y_pos = cur_node->position()[1] + ydisp;
      output_file << cur_x_pos << " " << cur_y_pos << " " << xdisp << " " << ydisp 
                  << cur_node->velocity()[0] << cur_node->velocity()[1]
                  << ((cur_node->material())->damageModel())->damageIndex() << std::endl;
    }
  }

  output_file.close();
 
  // Incremenet the output file count
  incrementOutputFileCount();
}

namespace Emu2DC {

  std::ostream& operator<<(std::ostream& out, const Output& output)
  {
    out.setf(std::ios::floatfield);
    out.precision(6);
    out << "Output dir = " << output.d_output_folder_name << " Output file = " << output.d_output_file_name
        << std::endl;
    out << "  Output iteration interval = " << output.d_output_iter_interval << std::endl;
    return out;
  }

}
