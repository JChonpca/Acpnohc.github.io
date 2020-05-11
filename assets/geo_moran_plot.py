# -*- coding: utf-8 -*-
"""
Created on Sat May  9 03:41:41 2020

@author: JChonpca_Huang
"""

from pyecharts import Map, Geo



def global_plot(province,values,name,year,min_num,max_num):
    
    #provice=list(province_distribution.keys())
    #values=list(province_distribution.values())
    
    # maptype='china' 只显示全国直辖市和省级
    # 数据只能是省名和直辖市的名称
    map = Map(name,year,width=1200, height=600)
    map.add('', province, values, visual_range=[0, 2],  maptype='china', is_visualmap=True,
        visual_text_color='#000')
    map.show_config()
    map.render(path=str(name) + '_' + str(year) +".html")
    
data = indexs[40:60,:]

province =[]

for i in sea_area:
    
    province.append(provin[i])
    
    

for plot in range(20):   #20X
    
    name_tmp = '1_'
    
    year_tmp = str(plot)
        
    values = data[plot,:].reshape(-1,1).tolist()
        
    global_plot(province,values,name_tmp,year_tmp,0,2)    
