package phoenixcenter.metaproteomics.entity;

import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.NoArgsConstructor;

import java.util.List;

@Data
@NoArgsConstructor
@AllArgsConstructor
@Builder
public class QuantPeptide {

    private String sequence;

    private List<Double> quantValList;
}
