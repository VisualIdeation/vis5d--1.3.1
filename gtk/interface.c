/*
 * DO NOT EDIT THIS FILE - it is generated by Glade.
 */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>

#include "support_cb.h"
#include "interface.h"
#include "support.h"

GtkWidget*
create_fileselection1 (void)
{
  GtkWidget *fileselection1;
  GtkWidget *ok_button1;
  GtkWidget *cancel_button1;

  fileselection1 = gtk_file_selection_new (_("Select File"));
  gtk_widget_set_name (fileselection1, "fileselection1");
  gtk_object_set_data (GTK_OBJECT (fileselection1), "fileselection1", fileselection1);
  gtk_container_set_border_width (GTK_CONTAINER (fileselection1), 10);
  gtk_window_set_modal (GTK_WINDOW (fileselection1), TRUE);

  ok_button1 = GTK_FILE_SELECTION (fileselection1)->ok_button;
  gtk_widget_set_name (ok_button1, "ok_button1");
  gtk_object_set_data (GTK_OBJECT (fileselection1), "ok_button1", ok_button1);
  gtk_widget_show (ok_button1);
  GTK_WIDGET_SET_FLAGS (ok_button1, GTK_CAN_DEFAULT);

  cancel_button1 = GTK_FILE_SELECTION (fileselection1)->cancel_button;
  gtk_widget_set_name (cancel_button1, "cancel_button1");
  gtk_object_set_data (GTK_OBJECT (fileselection1), "cancel_button1", cancel_button1);
  gtk_widget_show (cancel_button1);
  GTK_WIDGET_SET_FLAGS (cancel_button1, GTK_CAN_DEFAULT);

  gtk_signal_connect (GTK_OBJECT (ok_button1), "clicked",
                      GTK_SIGNAL_FUNC (on_fileselect_ok),
                      NULL);
  gtk_signal_connect (GTK_OBJECT (cancel_button1), "clicked",
                      GTK_SIGNAL_FUNC (on_fileselect_cancel),
                      NULL);

  return fileselection1;
}

GtkWidget*
create_fontselectiondialog1 (void)
{
  GtkWidget *fontselectiondialog1;
  GtkWidget *fontselectionokaybutton;
  GtkWidget *fontselectioncancelbutton;
  GtkWidget *fontselectionapplybutton;

  fontselectiondialog1 = gtk_font_selection_dialog_new (_("Select Font"));
  gtk_widget_set_name (fontselectiondialog1, "fontselectiondialog1");
  gtk_object_set_data (GTK_OBJECT (fontselectiondialog1), "fontselectiondialog1", fontselectiondialog1);
  gtk_container_set_border_width (GTK_CONTAINER (fontselectiondialog1), 4);
  gtk_window_set_policy (GTK_WINDOW (fontselectiondialog1), FALSE, TRUE, TRUE);

  fontselectionokaybutton = GTK_FONT_SELECTION_DIALOG (fontselectiondialog1)->ok_button;
  gtk_widget_set_name (fontselectionokaybutton, "fontselectionokaybutton");
  gtk_object_set_data (GTK_OBJECT (fontselectiondialog1), "fontselectionokaybutton", fontselectionokaybutton);
  gtk_widget_show (fontselectionokaybutton);
  GTK_WIDGET_SET_FLAGS (fontselectionokaybutton, GTK_CAN_DEFAULT);

  fontselectioncancelbutton = GTK_FONT_SELECTION_DIALOG (fontselectiondialog1)->cancel_button;
  gtk_widget_set_name (fontselectioncancelbutton, "fontselectioncancelbutton");
  gtk_object_set_data (GTK_OBJECT (fontselectiondialog1), "fontselectioncancelbutton", fontselectioncancelbutton);
  gtk_widget_show (fontselectioncancelbutton);
  GTK_WIDGET_SET_FLAGS (fontselectioncancelbutton, GTK_CAN_DEFAULT);

  fontselectionapplybutton = GTK_FONT_SELECTION_DIALOG (fontselectiondialog1)->apply_button;
  gtk_widget_set_name (fontselectionapplybutton, "fontselectionapplybutton");
  gtk_object_set_data (GTK_OBJECT (fontselectiondialog1), "fontselectionapplybutton", fontselectionapplybutton);
  gtk_widget_show (fontselectionapplybutton);
  GTK_WIDGET_SET_FLAGS (fontselectionapplybutton, GTK_CAN_DEFAULT);

  gtk_signal_connect (GTK_OBJECT (fontselectionokaybutton), "clicked",
                      GTK_SIGNAL_FUNC (on_fontselectionbutton_clicked),
                      0);
  gtk_signal_connect (GTK_OBJECT (fontselectioncancelbutton), "clicked",
                      GTK_SIGNAL_FUNC (on_fontselectionbutton_clicked),
                      2);
  gtk_signal_connect (GTK_OBJECT (fontselectionapplybutton), "clicked",
                      GTK_SIGNAL_FUNC (on_fontselectionbutton_clicked),
                      1);

  return fontselectiondialog1;
}

GtkWidget*
create_colorselectiondialog1 (void)
{
  GtkWidget *colorselectiondialog1;
  GtkWidget *ok_button2;
  GtkWidget *cancel_button2;
  GtkWidget *help_button1;

  colorselectiondialog1 = gtk_color_selection_dialog_new (_("Select Color"));
  gtk_widget_set_name (colorselectiondialog1, "colorselectiondialog1");
  gtk_object_set_data (GTK_OBJECT (colorselectiondialog1), "colorselectiondialog1", colorselectiondialog1);
  gtk_container_set_border_width (GTK_CONTAINER (colorselectiondialog1), 10);
  gtk_window_set_modal (GTK_WINDOW (colorselectiondialog1), TRUE);

  ok_button2 = GTK_COLOR_SELECTION_DIALOG (colorselectiondialog1)->ok_button;
  gtk_widget_set_name (ok_button2, "ok_button2");
  gtk_object_set_data (GTK_OBJECT (colorselectiondialog1), "ok_button2", ok_button2);
  gtk_widget_show (ok_button2);
  GTK_WIDGET_SET_FLAGS (ok_button2, GTK_CAN_DEFAULT);

  cancel_button2 = GTK_COLOR_SELECTION_DIALOG (colorselectiondialog1)->cancel_button;
  gtk_widget_set_name (cancel_button2, "cancel_button2");
  gtk_object_set_data (GTK_OBJECT (colorselectiondialog1), "cancel_button2", cancel_button2);
  gtk_widget_show (cancel_button2);
  GTK_WIDGET_SET_FLAGS (cancel_button2, GTK_CAN_DEFAULT);

  help_button1 = GTK_COLOR_SELECTION_DIALOG (colorselectiondialog1)->help_button;
  gtk_widget_set_name (help_button1, "help_button1");
  gtk_object_set_data (GTK_OBJECT (colorselectiondialog1), "help_button1", help_button1);
  gtk_widget_show (help_button1);
  GTK_WIDGET_SET_FLAGS (help_button1, GTK_CAN_DEFAULT);

  gtk_signal_connect (GTK_OBJECT (ok_button2), "clicked",
                      GTK_SIGNAL_FUNC (on_ColorSelectionOk_clicked),
                      NULL);
  gtk_signal_connect (GTK_OBJECT (cancel_button2), "clicked",
                      GTK_SIGNAL_FUNC (on_ColorSelectionCancel_clicked),
                      NULL);

  return colorselectiondialog1;
}

GtkWidget*
create_VerifyDialog (void)
{
  GtkWidget *VerifyDialog;
  GtkWidget *dialog_vbox1;
  GtkWidget *label1;
  GtkWidget *dialog_action_area1;
  GtkWidget *hbuttonbox1;
  GtkWidget *button1;
  GtkWidget *button2;

  VerifyDialog = gtk_dialog_new ();
  gtk_widget_set_name (VerifyDialog, "VerifyDialog");
  gtk_object_set_data (GTK_OBJECT (VerifyDialog), "VerifyDialog", VerifyDialog);
  gtk_window_set_title (GTK_WINDOW (VerifyDialog), _("Vis5d+"));
  gtk_window_set_policy (GTK_WINDOW (VerifyDialog), TRUE, TRUE, FALSE);

  dialog_vbox1 = GTK_DIALOG (VerifyDialog)->vbox;
  gtk_widget_set_name (dialog_vbox1, "dialog_vbox1");
  gtk_object_set_data (GTK_OBJECT (VerifyDialog), "dialog_vbox1", dialog_vbox1);
  gtk_widget_show (dialog_vbox1);

  label1 = gtk_label_new (_("Overwrite existing file ?"));
  gtk_widget_set_name (label1, "label1");
  gtk_widget_ref (label1);
  gtk_object_set_data_full (GTK_OBJECT (VerifyDialog), "label1", label1,
                            (GtkDestroyNotify) gtk_widget_unref);
  gtk_widget_show (label1);
  gtk_box_pack_start (GTK_BOX (dialog_vbox1), label1, FALSE, FALSE, 0);
  gtk_label_set_line_wrap (GTK_LABEL (label1), TRUE);
  gtk_misc_set_padding (GTK_MISC (label1), 5, 30);

  dialog_action_area1 = GTK_DIALOG (VerifyDialog)->action_area;
  gtk_widget_set_name (dialog_action_area1, "dialog_action_area1");
  gtk_object_set_data (GTK_OBJECT (VerifyDialog), "dialog_action_area1", dialog_action_area1);
  gtk_widget_show (dialog_action_area1);
  gtk_container_set_border_width (GTK_CONTAINER (dialog_action_area1), 10);

  hbuttonbox1 = gtk_hbutton_box_new ();
  gtk_widget_set_name (hbuttonbox1, "hbuttonbox1");
  gtk_widget_ref (hbuttonbox1);
  gtk_object_set_data_full (GTK_OBJECT (VerifyDialog), "hbuttonbox1", hbuttonbox1,
                            (GtkDestroyNotify) gtk_widget_unref);
  gtk_widget_show (hbuttonbox1);
  gtk_box_pack_start (GTK_BOX (dialog_action_area1), hbuttonbox1, TRUE, TRUE, 0);

  button1 = gtk_button_new_with_label (_("Okay"));
  gtk_widget_set_name (button1, "button1");
  gtk_widget_ref (button1);
  gtk_object_set_data_full (GTK_OBJECT (VerifyDialog), "button1", button1,
                            (GtkDestroyNotify) gtk_widget_unref);
  gtk_widget_show (button1);
  gtk_container_add (GTK_CONTAINER (hbuttonbox1), button1);
  GTK_WIDGET_SET_FLAGS (button1, GTK_CAN_DEFAULT);

  button2 = gtk_button_new_with_label (_("Cancel"));
  gtk_widget_set_name (button2, "button2");
  gtk_widget_ref (button2);
  gtk_object_set_data_full (GTK_OBJECT (VerifyDialog), "button2", button2,
                            (GtkDestroyNotify) gtk_widget_unref);
  gtk_widget_show (button2);
  gtk_container_add (GTK_CONTAINER (hbuttonbox1), button2);
  GTK_WIDGET_SET_FLAGS (button2, GTK_CAN_DEFAULT);

  return VerifyDialog;
}

