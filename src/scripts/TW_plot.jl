# 绘制Tianwen相关图像

module TW_plot
using Dates
using DelimitedFiles
using Interpolations
using LinearAlgebra
using DSP
using Wavelets
using CairoMakie
using GeometryBasics
using LaTeXStrings

export bowshock, magnetopause

const Rm = 3389.5 # 火星半径，单位km

"""
    bowshock(x)
火星弓激波模型
参数:
    x: 火星固连坐标系X坐标 (Rm)
返回值:
    弓激波半径 (Rm), Inf64 表示超出模型范围
"""
function bowshock(xshock)
    xF = 0.55
    ϵ = 1.05
    L = 2.10
    temp = (ϵ^2-1.0)*(xshock-xF)^2 - 2ϵ*L*(xshock-xF) + L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

"""
    magnetopause(x)
火星磁层顶模型
参数:
    x: 火星固连坐标系X坐标 (Rm)
返回值:
    磁层顶半径 (Rm), Inf64 表示超出模型范围
"""
function magnetopause(xmp)
    rSD = 1.33
    xF = 0.86
    ϵ = 0.92
    L = 0.90
    temp = (ϵ^2-1.0)*(xmp-xF)^2 - 2ϵ*L*(xmp-xF) + L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

end

