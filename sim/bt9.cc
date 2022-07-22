#include "bt9.h"

namespace bt9 {

/*!
     * \brief Overloaded output operator for BrBehavior::Direction class
     */
     std::ostream & operator<<(std::ostream & os, const BrBehavior::Direction & dir) 
    {
        os << StrEnumMap<BrBehavior::Direction>::convertEnumToStr(dir);
        return os;
    }

    /*!
     * \brief Overloaded output operator for BrBehavior::Indirectness class
     */
     std::ostream & operator<<(std::ostream & os, const BrBehavior::Indirectness & indir) 
    {
        os << StrEnumMap<BrBehavior::Indirectness>::convertEnumToStr(indir);
        return os;
    }

    /*!
     * \brief Overloaded output operator for BrBehavior class
     */
     std::ostream & operator<<(std::ostream & os, const BrBehavior & br_behav)
    {
        os << "behavior: " << std::setw(4) << br_behav.direction << "+" << br_behav.indirectness;

        return os;
    }

    /*!
     * \brief Overloaded output operator for BrClass::Type class
     */
     std::ostream & operator<<(std::ostream & os, const BrClass::Type & type)
    {
        os << StrEnumMap<BrClass::Type>::convertEnumToStr(type);
        return os;
    }

    /*!
     * \brief Overloaded output operator for BrClass::Directness class
     */
     std::ostream & operator<<(std::ostream & os, const BrClass::Directness & dir) 
    {
        os << StrEnumMap<BrClass::Directness>::convertEnumToStr(dir);
        return os;
    }
    
    /*!
     * \brief Overloaded output operator for BrClass::Conditionality class
     */
     std::ostream & operator<<(std::ostream & os, const BrClass::Conditionality & cond) 
    {
        os << StrEnumMap<BrClass::Conditionality>::convertEnumToStr(cond);
        return os;
    }

    /*!
     * \brief Overloaded output operator for BrClass class
     */
     std::ostream & operator<<(std::ostream & os, const BrClass & br_class) 
    {
        os << "class: ";
    
        os << std::setw(4) << br_class.type << "+" 
           << br_class.directness << "+" 
           << br_class.conditionality;

        return os;
    }

     /// Indicate if the static type (i.e. BrClass::Type) of this branch matches with the user-provided string
        bool BasicNodeRecord::brClassTypeIs(const std::string type_str) const
        {
            return (StrEnumMap<BrClass::Type>::convertEnumToStr(br_class_.type) == type_str);
        }

        /// Indicate if the static directness (i.e. BrClass::Directness) of this branch matches with the user-provided string
        bool BasicNodeRecord::brClassDirectnessIs(const std::string type_str) const
        {
            return (StrEnumMap<BrClass::Directness>::convertEnumToStr(br_class_.directness) == type_str);
        }

        /// Indicate if the static conditionality (i.e. BrClass::Conditionality) of this branch matches with the user-provided string
        bool BasicNodeRecord::brClassConditionalityIs(const std::string type_str) const
        {
            return (StrEnumMap<BrClass::Conditionality>::convertEnumToStr(br_class_.conditionality) == type_str);
        }

        /// Indicate if the dynamic direction (i.e. BrBehavior::Direction) of this branch matches with the user-provided string
        bool BasicNodeRecord::brBehaviorDirectionIs(const std::string type_str) const
        {
            return (StrEnumMap<BrBehavior::Direction>::convertEnumToStr(br_behavior_.direction) == type_str);
        }

        /// Indicate if the dynamic indirectness (i.e. BrBehavior::Indirectness) of this branch matches with the user-provided string
        bool BasicNodeRecord::brBehaviorIndirectnessIs(const std::string type_str) const
        {
            return (StrEnumMap<BrBehavior::Indirectness>::convertEnumToStr(br_behavior_.indirectness) == type_str);
        }




    /*!
     * \brief Specialize the lookup function of StrEnumMap template class with type 'BrBehavior::Indirectness'
     * \note The static lookup table defined in this function need to be updated every time
     *       changes are made on BrBehavior::Indirectness
     */
    template <>
    const EnumToStrMapType<BrBehavior::Indirectness> &StrEnumMap<BrBehavior::Indirectness>::getEnumToStrMap()
    {
        static const EnumToStrMapType<BrBehavior::Indirectness> map =
            {
                {BrBehavior::Indirectness::INDIRECT, "IND"},
                {BrBehavior::Indirectness::DIRECT, "DIR"}};

        return map;
    }

    /*!
     * \brief Specialize the lookup function of StrEnumMap template class with type 'BrClass::Type'
     * \note The static lookup table defined in this function need to be updated every time
     *       changes are made on BrClass::Type
     */
    template <>
    const EnumToStrMapType<BrClass::Type> &StrEnumMap<BrClass::Type>::getEnumToStrMap()
    {
        static const EnumToStrMapType<BrClass::Type> map =
            {
                {BrClass::Type::UNKNOWN, "N/A"},
                {BrClass::Type::RET, "RET"},
                {BrClass::Type::JMP, "JMP"},
                {BrClass::Type::CALL, "CALL"}};

        return map;
    }

    /*!
     * \brief Specialize the lookup function of StrEnumMap template class with type 'BrClass::Directness'
     * \note The static lookup table defined in this function need to be updated every time
     *       changes are made on BrClass::Directness
     */
    template <>
    const EnumToStrMapType<BrClass::Directness> &StrEnumMap<BrClass::Directness>::getEnumToStrMap()
    {
        static const EnumToStrMapType<BrClass::Directness> map =
            {
                {BrClass::Directness::UNKNOWN, "N/A"},
                {BrClass::Directness::DIRECT, "DIR"},
                {BrClass::Directness::INDIRECT, "IND"}};

        return map;
    }

    /*!
     * \brief Specialize the lookup function of StrEnumMap template class with type 'BrClass::Conditionality'
     * \note The static lookup table defined in this function need to be updated every time
     *       changes are made on BrClass::Conditionality
     */
    template <>
    const EnumToStrMapType<BrClass::Conditionality> &StrEnumMap<BrClass::Conditionality>::getEnumToStrMap()
    {
        static const EnumToStrMapType<BrClass::Conditionality> map =
            {
                {BrClass::Conditionality::UNKNOWN, "N/A"},
                {BrClass::Conditionality::CONDITIONAL, "CND"},
                {BrClass::Conditionality::UNCONDITIONAL, "UCD"}};

        return map;
    }
/*!
     * \brief Specialize the lookup function of StrEnumMap template class with type 'BrBehavior::Direction'
     * \note The static lookup table defined in this function need to be updated every time
     *       changes are made on BrBehavior::Direction
     */
    template <>
    const EnumToStrMapType<BrBehavior::Direction> &StrEnumMap<BrBehavior::Direction>::getEnumToStrMap()
    {
        static const EnumToStrMapType<BrBehavior::Direction> map =
            {
                {BrBehavior::Direction::AT, "AT"},
                {BrBehavior::Direction::ANT, "ANT"},
                {BrBehavior::Direction::DYN, "DYN"}};

        return map;
    }


    /*!
     * \brief Parse BrBehavior
     *
     * This is used by the BT9 reader library to parse the dynamic branch behavior
     * in the BT9 tracefile
     */
    void BrBehavior::parseBrBehavior(const std::string &str)
    {
        std::string token = str;
        std::replace(token.begin(), token.end(), '+', ' ');
        std::stringstream ss(token);

        while (ss >> token)
        {
            auto it1 = StrEnumMap<BrBehavior::Direction>::getStrToEnumMap().find(token);
            if (it1 != StrEnumMap<BrBehavior::Direction>::getStrToEnumMap().end())
            {
                direction = it1->second;
                continue;
            }

            auto it2 = StrEnumMap<BrBehavior::Indirectness>::getStrToEnumMap().find(token);
            if (it2 != StrEnumMap<BrBehavior::Indirectness>::getStrToEnumMap().end())
            {
                indirectness = it2->second;
                continue;
            }

            // Note: it1 and it2 are of different types.

            throw std::invalid_argument("Invalid token detected!");
        }
    }

    /*!
     * \brief Parse BrClass
     *
     * This is used by the BT9 reader library to parse the static branch categories
     * in the BT9 tracefile
     */
    void BrClass::parseBrClass(const std::string &str)
    {
        std::string token = str;
        std::replace(token.begin(), token.end(), '+', ' ');
        std::stringstream ss(token);

        while (ss >> token)
        {
            auto it1 = StrEnumMap<BrClass::Type>::getStrToEnumMap().find(token);
            if (it1 != StrEnumMap<BrClass::Type>::getStrToEnumMap().end())
            {
                type = it1->second;
                continue;
            }

            auto it2 = StrEnumMap<BrClass::Directness>::getStrToEnumMap().find(token);
            if (it2 != StrEnumMap<BrClass::Directness>::getStrToEnumMap().end())
            {
                directness = it2->second;
                continue;
            }

            auto it3 = StrEnumMap<BrClass::Conditionality>::getStrToEnumMap().find(token);
            if (it3 != StrEnumMap<BrClass::Conditionality>::getStrToEnumMap().end())
            {
                conditionality = it3->second;
                continue;
            }

            throw std::invalid_argument("Invalid token detected!");
        }

        // Some corner cases:
        if (type == BrClass::Type::RET)
        {
            directness = BrClass::Directness::INDIRECT;
        }
    }

    /*
     * \brief Convert a value of enum class type to std::string
     *
     * This is basically used by the overloaded output operator of enum class types
     */
    template <typename EnumType>
    const std::string &StrEnumMap<EnumType>::convertEnumToStr(const EnumType &enum_val)
    {
        auto it = StrEnumMap<EnumType>::getEnumToStrMap().find(enum_val);
        if (it != StrEnumMap<EnumType>::getEnumToStrMap().end())
        {
            return it->second;
        }
        else
        {
            throw std::invalid_argument("Undefined value detected!");
        }
    }
}
