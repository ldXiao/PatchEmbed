//
// Created by Lind Xiao on 7/4/19.
//
#pragma once
#ifndef OTMAPPING_EXPRESIONVALUE_H
#define OTMAPPING_EXPRESIONVALUE_H

#include <tinyexpr.h>
#include <nlohmann/json.hpp>
class ExpressionValue {
public:
    ~ExpressionValue();
    ExpressionValue();
    void init(const nlohmann::json &vals);
    void init(const double val);
    void init(const std::string &expr);

    double operator()(double x, double y) const;
    double operator()(double x, double y, double z) const;

private:
    struct Internal
    {
        double x, y, z;
    };

    te_expr *expr_;
    double value_;
    Internal *vals_;
};



#endif //OTMAPPING_EXPRESIONVALUE_H
