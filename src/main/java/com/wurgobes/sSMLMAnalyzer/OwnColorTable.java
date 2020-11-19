package com.wurgobes.sSMLMAnalyzer;

import java.awt.Color;
import java.io.IOException;

import net.imagej.lut.DefaultLUTService;
import net.imagej.lut.LUTService;

import net.imglib2.display.ColorTable;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * From https://forum.image.sc/t/plugins-java-how-to-get-a-color-from-current-specific-lut-using-an-index/4685/10
 * @author alex.vergara
 */
@Plugin(type = Command.class)
public class OwnColorTable{
    private Color[] CT;

    private final LUTService ls;

    public OwnColorTable(LUTService ls){
        this.ls = ls;
    }

    public void setLut(String colormap) throws IOException {
        ColorTable ct = ls.loadLUT(ls.findLUTs().get(colormap));
        int size = ct.getLength();
        CT = new Color[size];
        for (int i = 0; i < size; i++) {
            CT[i] = new Color(ct.get(ColorTable.RED, i), ct.get(ColorTable.GREEN, i), ct.get(ColorTable.BLUE, i));
        }
    }


    public Color getColor(float value, double min, double max){
        value = (float) (Math.max(Math.min((value - min)/(max - min), 1), 0) * 255); //remap from 0 to 1
        return CT[(int) value];
    }

    public String[] getLuts(){
        return ls.findLUTs().keySet().toArray(new String[0]);
    }

}