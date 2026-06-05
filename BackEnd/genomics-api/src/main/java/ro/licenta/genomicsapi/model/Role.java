package ro.licenta.genomicsapi.model;

/**
 * Role — rolurile disponibile în sistem.
 * USER: pacient/medic curent, vede doar propriile date
 * ADMIN: administrator, vede toate datele și utilizatorii
 */
public enum Role {
    USER,
    ADMIN
}