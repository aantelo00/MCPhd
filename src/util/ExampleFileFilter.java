package util;

/**********************************************
 * JACINTO D. CABRERA RODRIGUEZ               *
 *                Curso Java, Web y Applets   *
 * Basado en Jose F. Quesada & Jose A. Bernal *
 *  PRACTICA FINAL                            *
 * Fichero: ExampleFileFilter.java            *
 *********************************************/


import java.io.File;
import java.util.Hashtable;
import java.util.Enumeration;
import javax.swing.filechooser.*;

@SuppressWarnings("unchecked")
public class ExampleFileFilter extends FileFilter {

    @SuppressWarnings("rawtypes")
	private Hashtable filters = null;
    private String description = null;
    private String fullDescription = null;
    private boolean useExtensionsInDescription = true;
    @SuppressWarnings("unused")
	private String extension1;

    @SuppressWarnings("rawtypes")
	public ExampleFileFilter() {
	this.filters = new Hashtable();
    }

    public ExampleFileFilter(String extension) {
    	this(extension,null);
    	extension1 = extension;
    }

    public ExampleFileFilter(String extension, String description) {
	this();
	if(extension!=null) addExtension(extension);
 	if(description!=null) setDescription(description);
    }

    public ExampleFileFilter(String[] filters) {
	this(filters, null);
    }

    public ExampleFileFilter(String[] filters, String description) {
	this();
	for (int i = 0; i < filters.length; i++) {
	    addExtension(filters[i]);
	}
 	if(description!=null) setDescription(description);
    }

    public boolean accept(File f) {
	if(f != null) {
	    if(f.isDirectory()) {
		return true;
	    }
	    String extension = getExtension(f);
	    if(extension != null && filters.get(getExtension(f)) != null) {
		return true;
	    };
	}
	return false;
    }

    public String getExtension(File f) {
	if(f != null) {
	    String filename = f.getName();
	    int i = filename.lastIndexOf('.');
	    if(i>0 && i<filename.length()-1) {
		return filename.substring(i+1).toLowerCase();
	    };
	}
	return null;
    }

    @SuppressWarnings("rawtypes")
	public void addExtension(String extension) {
	if(filters == null) {
	    filters = new Hashtable(5);
	}
	filters.put(extension.toLowerCase(), this);
	fullDescription = null;
    }

    public String getDescription() {
	if(fullDescription == null) {
	    if(description == null || isExtensionListInDescription()) {
 		fullDescription = description==null ? "(" : description + " (";
		@SuppressWarnings("rawtypes")
		Enumeration extensions = filters.keys();
		if(extensions != null) {
		    fullDescription += "." + (String) extensions.nextElement();
		    while (extensions.hasMoreElements()) {
			fullDescription += ", " + (String) extensions.nextElement();
		    }
		}
		fullDescription += ")";
	    } else {
		fullDescription = description;
	    }
	}
	return fullDescription;
    }

    public void setDescription(String description) {
	this.description = description;
	fullDescription = null;
    }

    public void setExtensionListInDescription(boolean b) {
	useExtensionsInDescription = b;
	fullDescription = null;
    }

    public boolean isExtensionListInDescription() {
	return useExtensionsInDescription;
    }
}