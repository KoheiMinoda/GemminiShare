#include "mbconfig.h"           // This goes first in every *.c,*.cc file

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>

#include "dataman.h"
#include "userelem.h"
#include "strnode.h" // StructNodeを使用するために必要
#include "drive.h"   // DriveOwnerを使用するために必要
#include "module-catenary_lm.h" // モジュール登録用のヘッダ


// =================================================================
// === ModuleCatenaryLM クラス
// =================================================================

class ModuleCatenaryLM : virtual public Elem, public UserDefinedElem {
public:
    // コンストラクタとデストラクタ
    ModuleCatenaryLM(unsigned uLabel, const DofOwner *pDO, DataManager* pDM, MBDynParser& HP);
    virtual ~ModuleCatenaryLM(void);

    // MBDynインターフェース関数
    virtual void Output(OutputHandler& OH) const;
    virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

    SubVectorHandler&
    AssRes(SubVectorHandler& WorkVec,
        doublereal dCoef,
        const VectorHandler& XCurr, 
        const VectorHandler& XPrimeCurr);

    VariableSubMatrixHandler& 
    AssJac(VariableSubMatrixHandler& WorkMat,
        doublereal dCoef, 
        const VectorHandler& XCurr,
        const VectorHandler& XPrimeCurr);
    
    int iGetNumConnectedNodes(void) const;
    void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;

    // Initial系関数
    unsigned int iGetInitialNumDof(void) const;
    void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat, const VectorHandler& XCurr);
    SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

    // その他の必須仮想関数
    unsigned int iGetNumPrivData(void) const { return 0; };
    void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP, SimulationEntity::Hints *ph) { NO_OP; };
    std::ostream& Restart(std::ostream& out) const { return out << "# ModuleCatenaryLM: restart not implemented" << std::endl; };

private:
    std::vector<const StructNode*> m_nodes;

    // 物理パラメータ
    doublereal m_L_total;       // 係留索全長
    doublereal m_segment_length; // 1セグメントあたりの長さ
    doublereal m_EA;            // 軸剛性
    doublereal m_CA;            // 軸減衰
    doublereal m_rho_line;      // 水中線密度
    doublereal m_line_diameter; // 直径
    doublereal m_g_gravity;     // 重力加速度
    doublereal m_rho_water;     // 水の密度
    
    // 海底パラメータ
    doublereal m_seabed_z;      // 海底Z座標
    doublereal m_K_seabed;      // 海底剛性
    doublereal m_C_seabed;      // 海底減衰

    // 力のスケーリング用
    DriveOwner m_FSF;

    // 初期化フラグ
    bool m_bInitialStep;

    // 内部で使われる力計算関数
    Vec3 ComputeAxialForce(const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2) const;
    Vec3 ComputeGravityBuoyancy(doublereal segment_mass) const;
    Vec3 ComputeSeabedForce(const Vec3& position, const Vec3& velocity) const;

};

// =================================================================
// === コンストラクタの実装
// =================================================================

ModuleCatenaryLM::ModuleCatenaryLM(
    unsigned uLabel,
    const DofOwner *pDO,
    DataManager* pDM,
    MBDynParser& HP
)
: Elem(uLabel, flag(0)),
  UserDefinedElem(uLabel, pDO),
  m_bInitialStep(true)
{

    if (HP.IsKeyWord("help")) {}

    m_L_total = HP.GetReal();
    m_rho_line = HP.GetReal();
    m_EA = HP.GetReal();
    m_CA = HP.GetReal();
    m_line_diameter = HP.GetReal();
    m_g_gravity = HP.GetReal();
    m_rho_water = HP.GetReal();

    if (HP.IsKeyWord("seabed")) {
        m_seabed_z = HP.GetReal();
        m_K_seabed = HP.GetReal();
        m_C_seabed = HP.GetReal();
    } else {
        m_seabed_z = -320;
        m_K_seabed = 10000.0;
        m_C_seabed = 0.0;
    }

    int num_nodes = HP.GetInt();
    m_nodes.resize(num_nodes);

    for (int i = 0; i < num_nodes; ++i) {
        const Node* pNode = pDM->ReadNode(HP, Node::STRUCTURAL);
        
        if (!pNode) {
            silent_cerr("module-catenary_lm(" << GetLabel() << "): ERROR\n"
                << "    structural node expected but not found while reading node list.\n"
                << "    Check if all node labels in the .usr file are correctly defined in your .nod file.\n"
                << "    Error occurred reading " << i + 1 << "-th node in list at line: "
                << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        m_nodes[i] = dynamic_cast<const StructNode*>(pNode);

        if (!m_nodes[i]) {
            silent_cerr("module-catenary_lm(" << GetLabel() << "): ERROR\n"
                << "    Failed to cast node to StructNode.\n"
                << "    Node with label read from .usr file is not a structural node.\n"
                << "    Error occurred at line: " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    m_segment_length = m_L_total / (m_nodes.size() - 1);
    std::cout << "Calculated segment_length: " << m_segment_length << std::endl;

    // Force Scale Factor
    std::cout << "Setting up Force Scale Factor..." << std::endl;
    if (HP.IsKeyWord("force" "scale" "factor")) {
        m_FSF.Set(HP.GetDriveCaller());
    } else {
        m_FSF.Set(new OneDriveCaller);
    }

    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

    std::cout << "Final diagnosis:" << std::endl;
    int dynamic_count = 0, static_count = 0;
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (dynamic_cast<const DynamicStructNode*>(m_nodes[i])) {
            dynamic_count++;
        } else {
            static_count++;
        }
    }
    std::cout << "  Dynamic nodes: " << dynamic_count << std::endl;
    std::cout << "  Static nodes: " << static_count << std::endl;

    // ログファイルへの出力
    pDM->GetLogFile() << "catenary_lm: " << uLabel << " created, connected to " << m_nodes.size() << " nodes." << std::endl;
    pDM->GetLogFile() << "  Initial analysis enabled for static equilibrium calculation" << std::endl;
    pDM->GetLogFile() << "catenary_lm: " << uLabel << " parameters:" << std::endl;
    pDM->GetLogFile() << "  Total length: " << m_L_total << std::endl;
    pDM->GetLogFile() << "  Segment length: " << m_segment_length << std::endl;
    pDM->GetLogFile() << "  EA: " << m_EA << std::endl;
    pDM->GetLogFile() << "  CA: " << m_CA << std::endl;
    pDM->GetLogFile() << "  Seabed stiffness: " << m_K_seabed << std::endl;
}

// デストラクタ
ModuleCatenaryLM::~ModuleCatenaryLM(void) {}


// =================================================================
// === MBDynインターフェースの実装
// =================================================================

// 接続ノード数の報告
int ModuleCatenaryLM::iGetNumConnectedNodes(void) const {
    return m_nodes.size();
}

// 接続ノードのポインタリストをMBDynに渡す
void ModuleCatenaryLM::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
    connectedNodes.resize(m_nodes.size());
    for(size_t i = 0; i < m_nodes.size(); ++i) {
        connectedNodes[i] = m_nodes[i];
    }
}

// =================================================================
// === Initial系関数の実装
// =================================================================

// 初期自由度数の定義
unsigned int ModuleCatenaryLM::iGetInitialNumDof(void) const {
    int dynamic_node_count = 0;
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (dynamic_cast<const DynamicStructNode*>(m_nodes[i])) {
            dynamic_node_count++;
        }
    }
    return dynamic_node_count * 6;
}

// 初期作業領域のサイズ定義
void ModuleCatenaryLM::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
    int dynamic_node_count = 0;
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (dynamic_cast<const DynamicStructNode*>(m_nodes[i])) {
            dynamic_node_count++;
        }
    }
    *piNumRows = dynamic_node_count * 6;
    *piNumCols = dynamic_node_count * 6;
}

// 初期段階のヤコビアン行列
VariableSubMatrixHandler& ModuleCatenaryLM::InitialAssJac(
    VariableSubMatrixHandler& WorkMat, 
    const VectorHandler& XCurr) {
    
    int dynamic_node_count = 0;
    std::vector<bool> is_dynamic(m_nodes.size(), false);
    std::vector<int> node_to_dynamic_block(m_nodes.size(), -1);
    
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        const DynamicStructNode* dyn_node = dynamic_cast<const DynamicStructNode*>(m_nodes[i]);
        if (dyn_node != nullptr) {
            is_dynamic[i] = true;
            node_to_dynamic_block[i] = dynamic_node_count;
            dynamic_node_count++;
        }
    }

    if (dynamic_node_count == 0) {
        WorkMat.SetNullMatrix();
        return WorkMat;
    }

    FullSubMatrixHandler& WM = WorkMat.SetFull();
    WM.ResizeReset(dynamic_node_count * 6, dynamic_node_count * 6);

    for (size_t i=0; i<m_nodes.size(); ++i) {
        if (is_dynamic[i]) {
            integer iFirstIndex = m_nodes[i]->iGetFirstPositionIndex();
            int block = node_to_dynamic_block[i] * 6;
            
            for (int j=1; j<=6; j++) {
                WM.PutRowIndex(block + j, iFirstIndex + j);
                WM.PutColIndex(block + j, iFirstIndex + j);
            }
        }
    }

    doublereal dFSF = m_FSF.dGet();
    doublereal initial_EA = std::max(m_EA * 0.1, 1.0e6);

    std::cout << "ModuleCatenaryLM(" << GetLabel() << "): Initial Jacobian - using EA = " << initial_EA << std::endl;
    
    // セグメント間の初期剛性行列を構築
    for (size_t i=0; i<m_nodes.size() - 1; ++i) {
        
        const StructNode* node1 = m_nodes[i];
        const StructNode* node2 = m_nodes[i+1];
        
        const Vec3& x1 = node1->GetXCurr();
        const Vec3& x2 = node2->GetXCurr();
        Vec3 dx = x2 - x1;
        doublereal l_current = dx.Norm();
        
        if (l_current > 1e-12) {
            Vec3 t = dx / l_current;
            
            doublereal k_initial = initial_EA / l_current;
            Mat3x3 k_axial = Mat3x3(MatCrossCross, t, t) * k_initial * dFSF;

            doublereal k_geometric_scaler = initial_EA * 0.1 / (l_current * l_current);
            Mat3x3 k_geometric_mat = (Eye3 - Mat3x3(MatCrossCross, t, t)) * k_geometric_scaler * dFSF;
            Mat3x3 k_total = k_axial + k_geometric_mat;
            
            // 両方のノードが動的な場合
            if (is_dynamic[i] && is_dynamic[i+1]) {
                int block1 = node_to_dynamic_block[i] * 6 + 1;
                int block2 = node_to_dynamic_block[i+1] * 6 + 1;
                
                WM.Add(block1, block1, k_total);
                WM.Add(block2, block2, k_total);
                WM.Sub(block1, block2, k_total);
                WM.Sub(block2, block1, k_total);
                
            }
            //  フェアリーダーが静的な場合
            else if (is_dynamic[i] && !is_dynamic[i+1]) {
                int block1 = node_to_dynamic_block[i] * 6 + 1;
                WM.Add(block1, block1, k_total);
            } 
            // アンカーが静的な場合
            else if (!is_dynamic[i] && is_dynamic[i+1]) {
                int block2 = node_to_dynamic_block[i+1] * 6 + 1;
                WM.Add(block2, block2, k_total);
            }
        }
    }
    
    doublereal seabed_stiffness = std::max(m_K_seabed, 1.0e5);
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (!is_dynamic[i]) continue;
        
        const Vec3& pos = m_nodes[i]->GetXCurr();
        if (pos.dGet(3) < m_seabed_z && seabed_stiffness > 0.0) {
            int block = node_to_dynamic_block[i] * 3 + 3;
            WM.IncCoef(block, block, seabed_stiffness * dFSF);
        }
    }

    return WorkMat;
}

// 初期段階の残差ベクトル
SubVectorHandler& ModuleCatenaryLM::InitialAssRes(
    SubVectorHandler& WorkVec, 
    const VectorHandler& XCurr) {
    
    int dynamic_node_count = 0;
    std::vector<bool> is_dynamic(m_nodes.size(), false);
    std::vector<int> node_to_dynamic_block(m_nodes.size(), -1);
    
    // 動的ノードの特定
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        const DynamicStructNode* dyn_node = dynamic_cast<const DynamicStructNode*>(m_nodes[i]);
        if (dyn_node != nullptr) {
            is_dynamic[i] = true;
            node_to_dynamic_block[i] = dynamic_node_count;
            dynamic_node_count++;
        }
    }

    if (dynamic_node_count == 0) {
        WorkVec.ResizeReset(0);
        return WorkVec;
    }

    WorkVec.ResizeReset(dynamic_node_count * 6);

    // 動的ノードのインデックス設定
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (is_dynamic[i]) {
            integer iFirstIndex = m_nodes[i]->iGetFirstPositionIndex();
            int block = node_to_dynamic_block[i] * 6;
            
            for (int j = 1; j <= 6; j++) {
                WorkVec.PutRowIndex(block + j, iFirstIndex + j);
            }
        }
    }

    doublereal dFSF = m_FSF.dGet();
    double segment_mass = m_rho_line * m_segment_length;

    std::cout << "ModuleCatenaryLM(" << GetLabel() << "): Initial analysis - computing static equilibrium" << std::endl;

    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (!is_dynamic[i]) continue;
        
        const StructNode* current_node = m_nodes[i];
        const Vec3& x_current = current_node->GetXCurr();
        
        Vec3 F_total = Vec3(0.0, 0.0, 0.0);
        
        if (i > 0) {
            const StructNode* left_node = m_nodes[i-1];
            const Vec3& x_left = left_node->GetXCurr();
            
            Vec3 F_axial_left = ComputeAxialForce(x_left, x_current, Vec3(0,0,0), Vec3(0,0,0));
            F_total += F_axial_left;
        }
        
        if (i < m_nodes.size() - 1) {
            const StructNode* right_node = m_nodes[i+1];
            const Vec3& x_right = right_node->GetXCurr();
            
            Vec3 F_axial_right = ComputeAxialForce(x_current, x_right, Vec3(0,0,0), Vec3(0,0,0));
            F_total -= F_axial_right;
        }
        
        double segment_contribution = 0.0;
        if (i > 0) segment_contribution += 0.5;
        if (i < m_nodes.size() - 1) segment_contribution += 0.5;
        
        Vec3 F_gb = ComputeGravityBuoyancy(segment_mass);
        F_total += F_gb * segment_contribution;
        
        Vec3 F_seabed = ComputeSeabedForce(x_current, Vec3(0,0,0));
        F_total += F_seabed;
        
        int dof_start = node_to_dynamic_block[i] * 6 + 1;
        WorkVec.Add(dof_start, F_total * dFSF);
    }

    return WorkVec;
}

// =================================================================
// === 通常解析用関数の実装
// =================================================================

void ModuleCatenaryLM::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
    int dynamic_node_count = 0;

    for (size_t i=0; i<m_nodes.size(); ++i) {
        if (dynamic_cast<const DynamicStructNode*>(m_nodes[i])) {
            dynamic_node_count ++;
        }
    }

    *piNumRows = dynamic_node_count * 6;
    *piNumCols = dynamic_node_count * 6;
}

SubVectorHandler&
ModuleCatenaryLM::AssRes(
    SubVectorHandler& WorkVec, 
    doublereal dCoef,
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    if (m_bInitialStep) {
        std::cout << "ModuleCatenaryLM(" << GetLabel() << "): Transitioning from initial to dynamic analysis" << std::endl;
        m_bInitialStep = false;
    }

    int dynamic_node_count = 0;
    std::vector<bool> is_dynamic(m_nodes.size(), false);
    std::vector<int> node_to_dynamic_block(m_nodes.size(), -1);

    for (size_t i=0; i<m_nodes.size(); ++i) {
        const DynamicStructNode* dyn_node = dynamic_cast<const DynamicStructNode*>(m_nodes[i]);
        if (dyn_node != nullptr) {
            is_dynamic[i] = true;
            node_to_dynamic_block[i] = dynamic_node_count;
            dynamic_node_count ++;
        }
    }

    if (dynamic_node_count == 0) {
        WorkVec.ResizeReset(0);
        return WorkVec;
    }

    WorkVec.ResizeReset(dynamic_node_count*6);

    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (is_dynamic[i]) {
            integer iFirstIndex = m_nodes[i]->iGetFirstMomentumIndex();
            int dof_start = node_to_dynamic_block[i]*6;

            for (int j=1; j<=6; j++) {
                WorkVec.PutRowIndex(dof_start+j, iFirstIndex+j);
            }
        }
    }

    doublereal dFSF = m_FSF.dGet();
    double segment_mass = m_rho_line*m_segment_length;

    for (size_t i=0; i<m_nodes.size(); ++i) {
        if (!is_dynamic[i]) continue;

        const StructNode* current_node = m_nodes[i];
        const Vec3& x_current = current_node->GetXCurr();
        const Vec3& v_current = current_node->GetVCurr();

        Vec3 F_total = Vec3(0.0, 0.0, 0.0);

        if (i > 0) {
            const StructNode* left_node = m_nodes[i-1];
            const Vec3& x_left = left_node->GetXCurr();
            const Vec3& v_left = left_node->GetVCurr();
            
            Vec3 F_axial_left = ComputeAxialForce(x_left, x_current, v_left, v_current);
            F_total += F_axial_left;
        }

        if (i < m_nodes.size() - 1) {
            const StructNode* right_node = m_nodes[i+1];
            const Vec3& x_right = right_node->GetXCurr();
            const Vec3& v_right = right_node->GetVCurr();
            
            Vec3 F_axial_right = ComputeAxialForce(x_current, x_right, v_current, v_right);
            F_total -= F_axial_right;
        }

        double segment_contribution = 0.0;
        if (i > 0) segment_contribution += 0.5;
        if (i < m_nodes.size() - 1) segment_contribution += 0.5;
        
        Vec3 F_gb = ComputeGravityBuoyancy(segment_mass);
        F_total += F_gb * segment_contribution;
        
        Vec3 F_seabed = ComputeSeabedForce(x_current, v_current);
        F_total += F_seabed;
        
        int dof_start = node_to_dynamic_block[i] * 6 + 1;
        WorkVec.Add(dof_start, F_total * dFSF);

    }

    return WorkVec;
}

VariableSubMatrixHandler& 
ModuleCatenaryLM::AssJac(
    VariableSubMatrixHandler& WorkMat, 
    doublereal dCoef, 
    const VectorHandler& XCurr, 
    const VectorHandler& XPrimeCurr)
{
    int dynamic_node_count = 0;
    std::vector<bool> is_dynamic(m_nodes.size(), false);
    std::vector<int> node_to_dynamic_block(m_nodes.size(), -1);

    for (size_t i = 0; i < m_nodes.size(); ++i) {
        const DynamicStructNode* dyn_node = dynamic_cast<const DynamicStructNode*>(m_nodes[i]);
        if (dyn_node != nullptr) {
            is_dynamic[i] = true;
            node_to_dynamic_block[i] = dynamic_node_count;
            dynamic_node_count++;
        }
    }

    if (dynamic_node_count == 0) {
        WorkMat.SetNullMatrix();
        return WorkMat;
    }

    FullSubMatrixHandler& WM = WorkMat.SetFull();
    WM.ResizeReset(dynamic_node_count * 6, dynamic_node_count * 6);

    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (is_dynamic[i]) {
            integer iFirstIndex = m_nodes[i]->iGetFirstPositionIndex();
            int block = node_to_dynamic_block[i] * 6;
            
            for (int j = 1; j <= 6; j++) {
                WM.PutRowIndex(block + j, iFirstIndex + j);
                WM.PutColIndex(block + j, iFirstIndex + j);
            }
        }
    }

    doublereal dFSF = m_FSF.dGet();
    
    doublereal segment_mass = m_rho_line * m_segment_length;
    doublereal weight_per_segment = segment_mass * m_g_gravity;
    doublereal T_min = std::max(weight_per_segment, 10000.0);

    for (size_t i = 0; i < m_nodes.size() - 1; ++i) {

        const StructNode* node1 = m_nodes[i];
        const StructNode* node2 = m_nodes[i+1];

        const Vec3& x1 = node1->GetXCurr();
        const Vec3& x2 = node2->GetXCurr();
        Vec3 dx = x2 - x1;
        doublereal l_current = dx.Norm();
        
        if (l_current > 1e-12) {
            Vec3 t = dx / l_current;
            doublereal strain = (l_current - m_segment_length) / m_segment_length;
            Mat3x3 Total_stiffness;
            
            if (strain > 1e-6) {
                doublereal k_tangent = m_EA / l_current;
                Mat3x3 k_material = Mat3x3(MatCrossCross, t, t) * k_tangent;
                
                doublereal tension = std::max(m_EA * strain + T_min, T_min * 0.2);
                doublereal k_geometric_scaler = tension / (l_current * l_current);
                Mat3x3 k_geometric = (Eye3 - Mat3x3(MatCrossCross, t, t)) * k_geometric_scaler;

                Mat3x3 k_total = (k_material + k_geometric) * dFSF;
                Mat3x3 C_elem = k_material * (m_CA / m_EA) * dCoef * dFSF;
                Total_stiffness = k_total + C_elem;
                
            } else {
                doublereal k_tangent = T_min*0.1 / l_current;
                Total_stiffness = Mat3x3(MatCrossCross, t, t) * k_tangent * dFSF;
            }

            if (is_dynamic[i] && is_dynamic[i+1]) {
                int block1 = node_to_dynamic_block[i] * 6 + 1;
                int block2 = node_to_dynamic_block[i+1] * 6 + 1;
                
                WM.Add(block1, block1, Total_stiffness);
                WM.Add(block2, block2, Total_stiffness);
                WM.Sub(block1, block2, Total_stiffness);
                WM.Sub(block2, block1, Total_stiffness);
                
            } else if (is_dynamic[i] && !is_dynamic[i+1]) {
                int block1 = node_to_dynamic_block[i] * 6 + 1;
                WM.Add(block1, block1, Total_stiffness);
                
            } else if (!is_dynamic[i] && is_dynamic[i+1]) {
                int block2 = node_to_dynamic_block[i+1] * 6 + 1;
                WM.Add(block2, block2, Total_stiffness);
            }
        }
    }
    
    doublereal seabed_stiffness = std::min(m_K_seabed, 1.0e6);
    for (size_t i = 0; i < m_nodes.size(); ++i) {
        if (!is_dynamic[i]) continue;

        const Vec3& pos1 = m_nodes[i]->GetXCurr();
        if (pos1.dGet(3) < m_seabed_z && m_K_seabed > 0.0) {
            int block = node_to_dynamic_block[i] * 6 + 3; 
            WM.IncCoef(block, block, seabed_stiffness * dFSF);

            if (m_C_seabed > 0.0) {
                doublereal seabed_damping = std::min(m_C_seabed, 1.0e4) * dCoef;
                WM.IncCoef(block, block, seabed_damping * dFSF);
            }
        }
    }
        
    return WorkMat;
}

void ModuleCatenaryLM::Output(OutputHandler& OH) const
{
    if (bToBeOutput()) {
        if (OH.UseText(OutputHandler::LOADABLE)) {
            const Vec3& fairlead_pos = m_nodes.front()->GetXCurr();
            const Vec3& node1_pos = m_nodes[1]->GetXCurr();
            const Vec3& fairlead_vel = m_nodes.front()->GetVCurr();
            const Vec3& node1_vel = m_nodes[1]->GetVCurr();

            const Vec3& last_vnode_pos = m_nodes[m_nodes.size() - 2]->GetXCurr();
            const Vec3& anchor_pos = m_nodes.back()->GetXCurr();
            const Vec3& last_vnode_vel = m_nodes[m_nodes.size() - 2]->GetVCurr();
            const Vec3& anchor_vel = m_nodes.back()->GetVCurr();
            
            Vec3 F_tension_fairlead = ComputeAxialForce(fairlead_pos, node1_pos, fairlead_vel, node1_vel);
            Vec3 F_tension_anchor = ComputeAxialForce(last_vnode_pos, anchor_pos, last_vnode_vel, anchor_vel);

            OH.Loadable() << GetLabel()
                << " fairlead_tension " << F_tension_fairlead.Norm()
                << " anchor_tension " << F_tension_anchor.Norm()
                << std::endl;
        }
    }
}


// =================================================================
// === 内部の物理計算関数の実装
// =================================================================

Vec3 ModuleCatenaryLM::ComputeAxialForce(const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2) const {
    Vec3 dx = x2 - x1;
    doublereal l_current = dx.Norm();
    
    if (l_current < 1e-12) {
        return Vec3(0.0, 0.0, 0.0);
    }
    
    Vec3 t = dx / l_current;
    doublereal strain = (l_current - m_segment_length) / m_segment_length;

    doublereal segment_mass = m_rho_line * m_segment_length;
    doublereal weight_per_segment = segment_mass * m_g_gravity;
    doublereal T_min = std::max(weight_per_segment, 10000.0); 
    
    doublereal F_elastic = 0.0;
    if (strain > -0.05) { 
        F_elastic = m_EA * strain;
    } else {
        F_elastic = T_min;
    }
    
    Vec3 dv = v2 - v1;
    doublereal v_axial = dv.Dot(t);
    doublereal F_damping = m_CA * v_axial;
    
    doublereal F_total = F_elastic + F_damping;
    
    F_total = std::max(F_total, T_min);
    
    return t * F_total;
}

Vec3 ModuleCatenaryLM::ComputeGravityBuoyancy(doublereal segment_mass) const {
    doublereal weight = segment_mass * m_g_gravity;
    doublereal volume = (M_PI / 4.0) * m_line_diameter * m_line_diameter * m_segment_length;
    doublereal buoyancy = m_rho_water * volume * m_g_gravity;
    
    return Vec3(0.0, 0.0, -(weight - buoyancy));
}

Vec3 ModuleCatenaryLM::ComputeSeabedForce(const Vec3& position, const Vec3& velocity) const {
    Vec3 F_seabed(0.0, 0.0, 0.0);
    
    doublereal z = position.dGet(3);
    doublereal penetration = m_seabed_z - z;
    
    if (penetration > 0.0) {
        doublereal Fz = m_K_seabed * penetration - m_C_seabed * velocity.dGet(3);
        F_seabed = Vec3(0.0, 0.0, std::max(0.0, Fz)); // 上向きのみ
    }
    
    return F_seabed;
}


// =================================================================
// === モジュール登録のための定型句
// =================================================================

bool catenary_lm_set(void) {
    #ifdef DEBUG
        std::cerr << __FILE__ << ":" << __LINE__ << ":" << __PRETTY_FUNCTION__ << std::endl;
    #endif

    UserDefinedElemRead *rf = new UDERead<ModuleCatenaryLM>;

    if (!SetUDE("catenary_lm", rf)) {
        delete rf;
        return false;
    }

    return true;
}

#ifndef STATIC_MODULES
extern "C" {
    int module_init(const char *module_name, void *pdm, void *php) {
        if (!catenary_lm_set()) {
            silent_cerr("catenary_lm: "
                "module_init(" << module_name << ") "
                "failed" << std::endl);
            return -1;
        }
        return 0;
    }
}
#endif // ! STATIC_MODULES
