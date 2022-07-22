/*
 * Copyright 2015 Samsung Austin Semiconductor, LLC.
 */

/*!
 * \file    bt9.h
 * \brief   This is header file Branch Trace version 9 (BT9) related tools.
 */

#ifndef __BT9_H__
#define __BT9_H__

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <array>

namespace bt9
{
    using NodeTableHashKey = uint64_t;
    using EdgeTableHashKey = std::pair<uint64_t, uint64_t>;
}

namespace std
{
    template <>
    struct hash<bt9::EdgeTableHashKey>
    {
        size_t operator()(const bt9::EdgeTableHashKey &key) const
        {
            return hash<bt9::NodeTableHashKey>()(key.first) ^ hash<bt9::NodeTableHashKey>()(key.second);
        }
    };
}

namespace bt9
{
    /*!
     * \class BrBehavior
     * \brief Branch Behavior indicates dynamic observed attributes
     */
    class BrBehavior
    {
    public:
        enum class Direction
        {
            UNKNOWN,
            AT,  //!< Always Taken
            ANT, //!< Always Not-Taken
            DYN, //!< Dynamic
            // LOOP,
            // DAD,  //!< Dynamic Actually Direct
            // DAI,  //!< Dynamic Actually Indirect
            NUM_OF_BEHAV
        };

        enum class Indirectness
        {
            UNKNOWN,
            INDIRECT, //!< Behave Indirect
            DIRECT,   //!< Behave Direct
            NUM_OF_BEHAV
        };

        BrBehavior() = default;
        BrBehavior(const BrBehavior &) = default;

        void parseBrBehavior(const std::string &);

        Direction direction = BrBehavior::Direction::UNKNOWN;
        Indirectness indirectness = BrBehavior::Indirectness::UNKNOWN;
    };

    /*!
     * \class BrClass
     * \brief Branch Class indicates branch static attributes
     */
    class BrClass
    {
    public:
        enum class Type
        {
            UNKNOWN,
            RET,  //!< Return
            JMP,  //!< Jump
            CALL, //!< Call
        };

        enum class Directness
        {
            UNKNOWN,
            DIRECT,  //!< Direct
            INDIRECT //!< Indirect
        };

        enum class Conditionality
        {
            UNKNOWN,
            CONDITIONAL,  //!< Conditional
            UNCONDITIONAL //!< Unconditional
        };

        BrClass() = default;
        BrClass(const BrClass &) = default;

         void parseBrClass(const std::string &);

        Type type = BrClass::Type::UNKNOWN;
        Directness directness = BrClass::Directness::UNKNOWN;
        Conditionality conditionality = BrClass::Conditionality::UNKNOWN;
    };

    // Template alias C++-11 semantics
    template <typename T>
    using EnumToStrMapType = std::unordered_map<T, std::string>;

    template <typename T>
    using StrToEnumMapType = std::unordered_map<std::string, T>;

    /*!
     * \class StrEnumMap
     * \brief This is a helper class to perform conversion between user-defined enum class type and std::string
     */
    template <typename EnumType>
    class StrEnumMap
    {
    public:
        using EnumToStrMap = EnumToStrMapType<EnumType>;
        using StrToEnumMap = StrToEnumMapType<EnumType>;

        /*!
         * \brief Return the lookup table for enum class to string conversion
         *
         * \note This function needs to be specialized for every user-defined enum class type.
         *
         * \note A member or a member template of a class template may be explicitly specialized
         *       for a given implicit instantiation of the class template, even if the member or
         *       member template is defined in the class template definition.
         */
        static const EnumToStrMap &getEnumToStrMap();

        /*!
         * \brief Return the lookup table for string to enum class conversion
         */
        static const StrToEnumMap &getStrToEnumMap()
        {
            static const StrToEnumMap map = generateStrToEnumMap_();

            return map;
        }

        /*!
         * \brief Convert the given enum class instance to std::string
         */
        static const std::string &convertEnumToStr(const EnumType &);

    private:
        /*!
         * \brief This is the internal helper function to automatically generate the reverse lookup table
         *        from the enum to string map.
         */
        static StrToEnumMap generateStrToEnumMap_()
        {
            StrToEnumMap map;
            for (const auto &m : getEnumToStrMap())
            {
                map[m.second] = m.first;
            }

            return map;
        }
    };

    /*!
     * \brief Overloaded output operator for BrBehavior::Direction class
     */
    std::ostream &operator<<(std::ostream &os, const BrBehavior::Direction &dir);

    /*!
     * \brief Overloaded output operator for BrBehavior::Indirectness class
     */
    std::ostream &operator<<(std::ostream &os, const BrBehavior::Indirectness &indir);

    /*!
     * \brief Overloaded output operator for BrBehavior class
     */
    std::ostream &operator<<(std::ostream &os, const BrBehavior &br_behav);

    /*!
     * \brief Overloaded output operator for BrClass::Type class
     */
    std::ostream &operator<<(std::ostream &os, const BrClass::Type &type);

    /*!
     * \brief Overloaded output operator for BrClass::Directness class
     */
    std::ostream &operator<<(std::ostream &os, const BrClass::Directness &dir);

    /*!
     * \brief Overloaded output operator for BrClass::Conditionality class
     */
    std::ostream &operator<<(std::ostream &os, const BrClass::Conditionality &cond);

    /*!
     * \brief Overloaded output operator for BrClass class
     */
    std::ostream &operator<<(std::ostream &os, const BrClass &br_class);

    /////////////////////////////////////////////////////////////////////////////////////
    /// BT9 Header
    /////////////////////////////////////////////////////////////////////////////////////

    /*!
     * \class BasicHeader
     * \brief It stores the basic information of a BT9 header
     */
    class BasicHeader
    {
    public:
        /// BT9 file minor version number (starting from version 0)
        enum class BT9MinorVersionNum
        {
            VERSION_ORIGINAL = 0, //!< Original Version
            VERSION_MAX           //!< Invalid version number
        };

        BasicHeader() = default;
        BasicHeader(const BasicHeader &) = default;

        BasicHeader(BT9MinorVersionNum num,
                    bool has_phy_addr,
                    const std::string &md5,
                    const std::string &date,
                    std::string path) : version_num_(num),
                                        has_phy_addr_(has_phy_addr),
                                        md5sum_(md5),
                                        date_(date),
                                        original_tracefile_path_(path)
        {
        }

        const std::string &getOriginalTracefilePath() const { return original_tracefile_path_; }
        uint32_t getMinorVersionNum() const { return static_cast<uint32_t>(version_num_); }
        bool getHasPhyAddr() const { return has_phy_addr_; }
        const std::string &getMd5CheckSum() const { return md5sum_; }
        const std::string &getDate() const { return date_; }

        /// Extract trace name from original tracefile full path
        std::string extractTraceName() const
        {
            // Get Trace Name
            // Assume <file>.xxx.gz as input trace file name
            int lastindex = original_tracefile_path_.find_last_of(".");
            std::string trace_name = original_tracefile_path_.substr(0, lastindex);

            // Assume <file>.xxx after above string processing
            lastindex = trace_name.find_last_of(".");
            trace_name = trace_name.substr(0, lastindex);

            lastindex = trace_name.find_last_of("/\\");
            trace_name = trace_name.substr(lastindex + 1, -1);

            return trace_name;
        }

        virtual void printBasicInfo(std::ostream &os) const
        {
            os << "bt9_minor_version: " << getMinorVersionNum() << '\n'
               << "has_physical_address: " << getHasPhyAddr() << '\n'
               << "md5_checksum: " << getMd5CheckSum() << '\n'
               << "conversion_date: " << getDate() << '\n'
               << "original_stf_input_file: " << getOriginalTracefilePath() << '\n';
        }

        virtual ~BasicHeader() {}

    protected:
        /// BT9 file minor version
        BT9MinorVersionNum version_num_ = BT9MinorVersionNum::VERSION_ORIGINAL;

        /// Indicate if trace contains physical address info
        bool has_phy_addr_ = false;

        /// MD5 checksum of original trace file
        std::string md5sum_;

        /// Date and time of conversion
        std::string date_;

        /// Original tracefile
        std::string original_tracefile_path_;
    };

    /////////////////////////////////////////////////////////////////////////////////////
    /// Node Record
    /////////////////////////////////////////////////////////////////////////////////////

    /*!
     * \class BasicNodeRecord
     * \brief It stores the basic static encoded and dynamic observed attributes for each single branch
     */
    class BasicNodeRecord
    {
    public:
        BasicNodeRecord() = default;
        BasicNodeRecord(const BasicNodeRecord &) = default;

        BasicNodeRecord(uint32_t id,
                        uint64_t br_virtual_addr,
                        bool br_phy_addr_valid,
                        uint64_t br_phy_addr,
                        uint32_t opcode,
                        uint64_t opcode_size,
                        BrClass br_class,
                        BrBehavior br_behavior,
                        std::string mnemonic,
                        uint32_t br_taken_cnt,
                        uint32_t br_untaken_cnt) : id_(id),
                                                   br_virtual_addr_(br_virtual_addr),
                                                   br_phy_addr_valid_(br_phy_addr_valid),
                                                   br_phy_addr_(br_phy_addr),
                                                   opcode_(opcode),
                                                   opcode_size_(opcode_size),
                                                   br_class_(br_class),
                                                   br_behavior_(br_behavior),
                                                   mnemonic_(mnemonic),
                                                   br_taken_cnt_(br_taken_cnt),
                                                   br_untaken_cnt_(br_untaken_cnt)
        {
        }

        /// Get branch index
        uint32_t brNodeIndex() const { return id_; }

        /// Get branch virtual address
        uint64_t brVirtualAddr() const { return br_virtual_addr_; }

        /// Indicate if branch physical address is valid
        bool brPhyAddrIsValid() const { return br_phy_addr_valid_; }

        /// Get branch physical address (if it is valid)
        uint64_t brPhyAddr() const { return br_phy_addr_; }

        /// Get branch opcode
        uint32_t brOpcode() const { return opcode_; }

        /// Get branch opcode size
        uint64_t brOpcodeSize() const { return opcode_size_; }

        /// Get branch class (static encoding information)
        BrClass brClass() const { return br_class_; }

        /// Get branch behavior (dynamic execution information)
        BrBehavior brBehavior() const { return br_behavior_; }

        /// Get dynamically observed taken counts for this branch
        uint32_t brObservedTakenCnt() const { return br_taken_cnt_; }

        /// Get dynamically observed not-taken counts for this branch
        uint32_t brObservedNotTakenCnt() const { return br_untaken_cnt_; }

        /// Indicate if the static type (i.e. BrClass::Type) of this branch matches with the user-provided string
        bool brClassTypeIs(const std::string type_str) const;

        /// Indicate if the static directness (i.e. BrClass::Directness) of this branch matches with the user-provided string
        bool brClassDirectnessIs(const std::string type_str) const;

        /// Indicate if the static conditionality (i.e. BrClass::Conditionality) of this branch matches with the user-provided string
        bool brClassConditionalityIs(const std::string type_str) const;

        /// Indicate if the dynamic direction (i.e. BrBehavior::Direction) of this branch matches with the user-provided string
        bool brBehaviorDirectionIs(const std::string type_str) const;

        /// Indicate if the dynamic indirectness (i.e. BrBehavior::Indirectness) of this branch matches with the user-provided string
        bool brBehaviorIndirectnessIs(const std::string type_str) const;

        /*!
         * \brief print fixed fields of branch node record
         *
         * This is expected to be used for its derived class also
         */
        void printFixedFields(std::ostream &os) const
        {
            std::ios_base::fmtflags saveflags = os.flags();
            std::streamsize prec = os.precision();
            std::streamsize width = os.width();

            os << std::setw(4) << id_ << std::setw(1) << " "
               << std::hex << std::showbase << std::setw(20) << br_virtual_addr_
               << std::setw(1) << " ";

            if (br_phy_addr_valid_)
            {
                os << std::hex << std::showbase << std::setw(20) << br_phy_addr_
                   << std::setw(1) << " ";
            }
            else
            {
                os << std::setw(20) << "-" << std::setw(1) << " ";
            }

            os << std::hex << std::showbase << std::setw(16) << opcode_ << std::setw(1) << " "
               << std::dec << std::noshowbase << std::setw(4) << opcode_size_
               << std::setw(2) << "  ";

            os.flags(saveflags);
            os.precision(prec);
            os.width(width);
        }

        /*!
         * \brief Print optional fields of branch node record
         *
         * This could be overridden by its derived class
         */
        virtual void printOptionalFields(std::ostream &os) const
        {
            std::ios_base::fmtflags saveflags = os.flags();
            std::streamsize prec = os.precision();
            std::streamsize width = os.width();

            if (opcode_size_ == 0)
            {
                return;
            }

            // Print pre-defined key-value pair fields
            os << br_class_ << std::setw(2) << " "
               << br_behavior_ << std::setw(2) << " "
               << "taken_cnt: " << std::setw(8) << std::dec << br_taken_cnt_
               << std::setw(2) << " "
               << "not_taken_cnt: " << std::setw(8) << std::dec << br_untaken_cnt_;

            os.flags(saveflags);
            os.precision(prec);
            os.width(width);
        }

        virtual void printComments(std::ostream &os) const
        {
            (void)os;
        }

        virtual ~BasicNodeRecord() {}

        uint32_t edgeOut(uint32_t node_id, bool taken = true);

    protected:
        // Number of fixed value fields in node record
        static const uint8_t NUM_VALUE_FIELD_ = 5;

        // Pre-defined Fixed Fields
        uint32_t id_ = 0;
        uint64_t br_virtual_addr_ = 0;
        bool br_phy_addr_valid_ = 0;
        uint64_t br_phy_addr_ = 0;
        uint32_t opcode_ = 0;
        uint64_t opcode_size_ = 0;

        // Pre-defined Key-Value Pairs
        BrClass br_class_;
        BrBehavior br_behavior_;
        std::string mnemonic_;
        uint32_t br_taken_cnt_ = 0;
        uint32_t br_untaken_cnt_ = 0;
    };

    /*!
     * \brief Overloaded output operator for branch node record
     */
    inline std::ostream &operator<<(std::ostream &os, const BasicNodeRecord &rec)
    {
        rec.printFixedFields(os);
        os << "  ";
        rec.printOptionalFields(os);
        os << "  ";
        rec.printComments(os);

        return os;
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Edge Record
    /////////////////////////////////////////////////////////////////////////////////////

    /*!
     * \class BasicEdgeRecord
     * \brief It stores dynamically attributes for each edge (branch pairs) of control flow graph
     */
    class BasicEdgeRecord
    {
    public:
        BasicEdgeRecord() = default;
        BasicEdgeRecord(uint32_t id,
                        uint32_t src_node_id,
                        uint32_t dest_node_id,
                        bool is_taken_path,
                        uint64_t br_virtual_tgt,
                        bool br_phy_tgt_valid,
                        uint64_t br_phy_tgt,
                        uint64_t inst_cnt,
                        uint64_t observed_traverse_cnt) : id_(id),
                                                          src_node_id_(src_node_id),
                                                          dest_node_id_(dest_node_id),
                                                          is_taken_path_(is_taken_path),
                                                          br_virtual_tgt_(br_virtual_tgt),
                                                          br_phy_tgt_valid_(br_phy_tgt_valid),
                                                          br_phy_tgt_(br_phy_tgt),
                                                          inst_cnt_(inst_cnt),
                                                          observed_traverse_cnt_(observed_traverse_cnt)
        {
        }

        BasicEdgeRecord(const BasicEdgeRecord &) = default;

        /// Get edge index
        uint32_t edgeIndex() const { return id_; }

        /// Get source node index
        uint32_t srcNodeIndex() const { return src_node_id_; }

        /// Get destionation node index
        uint32_t destNodeIndex() const { return dest_node_id_; }

        /// Indicate if it is a taken path
        bool isTakenPath() const { return is_taken_path_; }

        /// Get branch target virtual address
        uint64_t brVirtualTarget() const { return br_virtual_tgt_; }

        /// Indicate if branch target physical address is valid
        bool brPhyTargetIsValid() const { return br_phy_tgt_valid_; }

        /// Get branch target physical address
        uint64_t brPhyTarget() const { return br_phy_tgt_; }

        /// Get non-branch instruction count on this path
        uint64_t nonBrInstCnt() const { return inst_cnt_; }

        /// Get the observed traverse count on this path
        uint64_t observedTraverseCnt() const { return observed_traverse_cnt_; }

        /*!
         * \brief print fixed fields of branch edge record
         *
         * This is expected to be used for its derived class also
         */
        void printFixedFields(std::ostream &os) const
        {
            std::ios_base::fmtflags saveflags = os.flags();
            std::streamsize prec = os.precision();
            std::streamsize width = os.width();

            os << std::setw(6) << id_ << std::setw(1) << " "
               << std::setw(6) << src_node_id_ << std::setw(1) << " "
               << std::setw(6) << dest_node_id_ << std::setw(1) << " ";

            if (is_taken_path_)
            {
                os << std::setw(8) << "T" << std::setw(1) << " ";
            }
            else
            {
                os << std::setw(8) << "N" << std::setw(1) << " ";
            }

            os << std::setw(20) << std::hex << std::showbase << br_virtual_tgt_
               << std::setw(1) << " ";

            if (br_phy_tgt_valid_)
            {
                os << std::setw(20) << std::hex << std::showbase << br_phy_tgt_
                   << std::setw(1) << " ";
            }
            else
            {
                os << std::setw(20) << "-" << std::setw(1) << " ";
            }

            os << std::setw(8) << std::dec << inst_cnt_ << std::setw(1) << " ";

            os.flags(saveflags);
            os.precision(prec);
            os.width(width);
        }

        /*!
         * \brief Print optional fields of branch edge record
         *
         * This could be overridden by its derived class
         */
        virtual void printOptionalFields(std::ostream &os) const
        {
            std::ios_base::fmtflags saveflags = os.flags();
            std::streamsize prec = os.precision();
            std::streamsize width = os.width();

            os << "traverse_cnt: " << std::setw(8) << std::dec
               << observed_traverse_cnt_ << std::setw(2) << " ";

            os.flags(saveflags);
            os.precision(prec);
            os.width(width);
        }

        virtual ~BasicEdgeRecord() {}

    protected:
        // Number of fixed fields in edge record
        static const uint8_t NUM_VALUE_FIELD_ = 7;

        // Fixed Fields
        uint32_t id_ = 0;
        uint32_t src_node_id_ = 0;
        uint32_t dest_node_id_ = 0;
        bool is_taken_path_ = 0;
        uint64_t br_virtual_tgt_ = 0;
        bool br_phy_tgt_valid_ = 0;
        uint64_t br_phy_tgt_ = 0;
        uint64_t inst_cnt_ = 0;

        // Key-Value Pairs
        uint64_t observed_traverse_cnt_ = 0;
    };

    /*!
     * \brief Overloaded output operator for branch edge record
     */
    inline std::ostream &operator<<(std::ostream &os, const BasicEdgeRecord &rec)
    {
        rec.printFixedFields(os);
        os << '\t';
        rec.printOptionalFields(os);

        return os;
    }



}

namespace std
{

    template <>
    struct hash<bt9::BrBehavior::Direction>
    {
        size_t operator()(const bt9::BrBehavior::Direction &key) const
        {
            return hash<uint8_t>()(static_cast<uint8_t>(key));
        }
    };

    template <>
    struct hash<bt9::BrBehavior::Indirectness>
    {
        size_t operator()(const bt9::BrBehavior::Indirectness &key) const
        {
            return hash<uint8_t>()(static_cast<uint8_t>(key));
        }
    };

    template <>
    struct hash<bt9::BrClass::Type>
    {
        size_t operator()(const bt9::BrClass::Type &key) const
        {
            return hash<uint8_t>()(static_cast<uint8_t>(key));
        }
    };

    template <>
    struct hash<bt9::BrClass::Directness>
    {
        size_t operator()(const bt9::BrClass::Directness &key) const
        {
            return hash<uint8_t>()(static_cast<uint8_t>(key));
        }
    };

    template <>
    struct hash<bt9::BrClass::Conditionality>
    {
        size_t operator()(const bt9::BrClass::Conditionality &key) const
        {
            return hash<uint8_t>()(static_cast<uint8_t>(key));
        }
    };

}

// __BT9_H__
#endif
