package com.wurgobes.sSMLMAnalyzer;

import java.io.File;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.Objects;

//https://stackabuse.com/example-loading-a-java-class-at-runtime/

public class ExtensionLoader<C> {

    public C LoadClass(String directory, String classpath, Class<C> parentClass) throws ClassNotFoundException {
        File pluginsDir = new File(directory);

        for (File jar : Objects.requireNonNull(pluginsDir.listFiles())) {
            try {
                ClassLoader loader = URLClassLoader.newInstance(
                        new URL[] { jar.toURI().toURL()},
                        getClass().getClassLoader()
                );
                Class<?> clazz = Class.forName(classpath, true, loader);

                Class<? extends C> newClass = clazz.asSubclass(parentClass);

                Constructor<? extends C> constructor = newClass.getConstructor();
                return constructor.newInstance();

            } catch (ClassNotFoundException e) {
                // There might be multiple JARs in the directory,
                // so keep looking
            } catch (MalformedURLException | NoSuchMethodException | InvocationTargetException | IllegalAccessException | InstantiationException e) {
                e.printStackTrace();
            }
        }
        throw new ClassNotFoundException("Class " + classpath
                + " wasn't found in directory " + System.getProperty("user.dir") + directory);
    }
}