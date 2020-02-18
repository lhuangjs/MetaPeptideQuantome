package phoenixcenter.metaproteomics.entity;

import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.NoArgsConstructor;

@Data
@NoArgsConstructor
@AllArgsConstructor
@Builder
public class UnipeptTaxon {

    private Integer id;

    private String name;

    private String rank;
}