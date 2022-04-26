package com.wurgobes.sSMLMAnalyzer;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import java.awt.Checkbox;
import java.awt.TextField;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.Vector;

public class GenericDialogPlusPlus extends GenericDialogPlus {
    public Button dummy = new Button();

    @Override
    public void actionPerformed(ActionEvent event) {
        super.actionPerformed(checkSettings(event));
    }

    public String validInputs() {
        Vector<Checkbox> checkboxes = super.getCheckboxes();
        Vector<TextField> strings = super.getStringFields();

        String filePath = strings.get(0).getText();
        String csv_target_dir = strings.get(1).getText();

        boolean saveSCV = checkboxes.get(0).getState();
        boolean visualisation = checkboxes.get(9).getState();
        boolean visualiseZOLA = checkboxes.get(10).getState();
        // Require input CSV
        if(filePath.equals("")){
            return"No input CSV was set";
        }

        // Require output directory if you want to save
        if(saveSCV && csv_target_dir.equals("")){
            return "Set saving to CSV but no filepath was provided.";
        }

        // Require user to either save or show results
        if(!(visualisation || saveSCV || visualiseZOLA)) {
            return"No output method of any sorts is selected.\nSelect either Visualisation or Save to CSV.";
        }

        return "OK";
    }

    private ActionEvent checkSettings(ActionEvent event) {
        if(event.getActionCommand().contains("OK")){
            String answer = validInputs();
            if(!answer.equals("OK")) {
                IJ.showMessage(answer);
                event.setSource(dummy);
            }
        }
        return event;
    }

    public GenericDialogPlusPlus(String title) {
        super(title);
    }
}
